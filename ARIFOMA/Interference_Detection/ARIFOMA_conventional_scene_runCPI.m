function outS = ARIFOMA_conventional_scene_runCPI( ...
        ax, numCPI, cfg, cIndex, ctx, V, I, ...
        fv_base, fi_base, schemeInterferer)
% ARIFOMA_CONVENTIONAL_SCENE_RUNCPI
%   Execute one CPI for a single victim corner using the
%   conventional (non-adaptive) FMCW interference analysis pipeline.
%
% PURPOSE
% -------
%   This function implements the core per-CPI analysis loop for the
%   CONVENTIONAL scheme and serves as the reference baseline against which
%   all mitigation strategies are compared.  For each CPI micro-step it:
%
%     1. Extracts victim and interferer waveform parameters.
%     2. Generates the victim LO signal (with waveform caching for static
%        schemes to avoid redundant regeneration).
%     3. Generates all interferer LO signals (also cached for static schemes).
%     4. Runs the analytic interference detector for each victim–interferer
%        pair and accumulates per-chirp and per-path-type metrics.
%     5. Packs all FOMs into the output struct outS.
%     6. Optionally renders the time–frequency spectrogram overlay on the
%        corner tile axes.
%
% WAVEFORM CACHING
% ----------------
%   LO signals are expensive to regenerate each CPI when the scheme is
%   static (carrier frequency does not change between CPIs).
%   ARIFOMA_scheme_is_dynamic determines whether a scheme requires fresh
%   generation each call:
%     Static  (cached):   CONVENTIONAL, COMPASS, FIXED PARAMETER POOL
%     Dynamic (no cache): FAR, RANDOM_FH_PER_CPI, RANDOM_FH_PER_CPI_INTERFERED,
%                         BATSINSPIRED_FH_PER_CPI
%
% INPUTS
%   ax                (axes)   — target corner tile axes handle
%   numCPI            (int)    — current CPI micro-step index (1-based);
%                                used in the tile title
%   cfg               (struct) — global configuration struct; required fields:
%                                  cfg.randstartTime           [s]
%                                  cfg.BADC_Hz                 [Hz]
%                                  cfg.boolIgnoreFlybackInterferer (logical)
%                                  cfg.CPI_duration_s          [s]
%                                  cfg.plotScenarioEnd         (logical)
%   cIndex            (int or string) — corner identifier for tile title
%   ctx               (struct) — current corner context; modified fields:
%                                  ctx.waveformCache — LO signal cache
%   V                 (struct) — victim sensor parameter row; required fields:
%                                  SensorID, Scheme, Waveform, B, T, Nch,
%                                  ChirpsPerStep, FreqIncrHz, Dir,
%                                  bandMin, bandMax,
%                                  T_flyback_default, T_flyback_stepped,
%                                  T_returnfly, Nseg, PermPattern
%   I                 (struct array, 1×nI) — interferer parameter rows;
%                                  same fields as V plus:
%                                  min_d1_m — nearest first-leg distance [m]
%                                  min_d2_m — nearest second-leg distance [m]
%                                  (NaN or 0 for LOS paths)
%   fv_base           (double) — victim carrier frequency F0 [Hz]
%   fi_base           (double, nI×1) — interferer carrier frequencies F0 [Hz]
%   schemeInterferer  (string) — active interferer scheme identifier;
%                                "NOTHING" suppresses interferer rendering
%
% OUTPUT
%   outS  (struct) — per-CPI result struct with fields:
%
%     Context
%       outS.ctx                          updated corner context (with cache)
%
%     CPI-level FOMs
%       outS.P_coll_per_n                 per-CPI interference ratio [0,1]
%       outS.total_interfered_chirps      number of victim chirps with any overlap
%
%     Normalised interference durations (T_I / T_chirp)
%       outS.all_interference_durations       combined LOS + REFL
%       outS.all_interference_durations_LOS   LOS only
%       outS.all_interference_durations_REFL  reflected only
%
%     Path distances
%       outS.all_interference_dist_LOS        LOS path distances [m]
%       outS.all_interference_dist_REFL       total reflected distances [m]
%       outS.all_interference_dist_REFL_d1    first reflected leg [m]
%       outS.all_interference_dist_REFL_d2    second reflected leg [m]
%
%     Per-chirp interferer counts
%       outS.all_nIntf_per_chirp_total    total interferers per chirp
%       outS.all_nIntf_per_chirp_LOS      LOS interferers per chirp
%       outS.all_nIntf_per_chirp_REFL     reflected interferers per chirp
%
%     BATS placeholder
%       outS.all_r_hat_sub                NaN (BATS not active in this strategy)
%
% NOTES
% -----
%   - A random start-time offset (uniform over [0, cfg.randstartTime]) is
%     applied to the victim's time grid to model unsynchronised radar
%     transmissions (concurrent transmission probability P_CT < 1).
%   - The collision mask collMask(k,j) is true if victim chirp k overlaps
%     with at least one interval from interferer j.  The per-CPI
%     interference ratio counts chirps where any(collMask(k,:)) is true.
%   - Visualisation is rendered only when cfg.plotScenarioEnd = true.
%     LOS interferers are shown in red; reflected interferers in grey.
%     Interference regions are shaded blue.
%
% SEE ALSO
%   ARIFOMA_analytic_detector,
%   ARIFOMA_generate_LO_Signal_by_Scheme,
%   ARIFOMA_scheme_is_dynamic,
%   ARIFOMA_run_strategy

% =========================================================================

%% ── Section A: Victim Waveform Parameters ────────────────────────────────
%   Unpack all victim waveform fields into local variables for readability
%   and to avoid repeated struct field lookups in the generation call.

Bv                  = V.B;
Tv                  = V.T;
Nch_V               = V.Nch;
chirpsPerStep_V     = V.ChirpsPerStep;
freqIncr_V          = V.FreqIncrHz;
dir_V               = V.Dir;
bandMin_V           = V.bandMin;
bandMax_V           = V.bandMax;
T_flyback_default_V = V.T_flyback_default;
T_flyback_stepped_V = V.T_flyback_stepped;
T_returnfly_V       = V.T_returnfly;
Nseg_V              = V.Nseg;
PermPattern_V       = V.PermPattern;

%% ── Section B: Interferer Path Distance Metadata ─────────────────────────
%   For each interferer, determine whether the dominant path to the victim
%   is LOS or reflected and store the relevant distances for FOM binning.

nI       = numel(I);
intf_meta = repmat(struct( ...
    'd_LOS',NaN, 'd_REFL_total',NaN, ...
    'd_REFL_d1',NaN, 'd_REFL_d2',NaN), 1, nI);

othersIDs_V = zeros(1, nI);

for m = 1:nI

    d1 = double(I(m).min_d1_m);
    d2 = double(I(m).min_d2_m);

    if isfinite(d2) && d2 > 0
        % ── Reflected path: both legs defined ────────────────────────────
        intf_meta(m).d_LOS        = NaN;
        intf_meta(m).d_REFL_d1    = d1;
        intf_meta(m).d_REFL_d2    = d2;
        intf_meta(m).d_REFL_total = d1 + d2;
    else
        % ── LOS path: only first leg defined ─────────────────────────────
        intf_meta(m).d_LOS        = d1;
        intf_meta(m).d_REFL_d1    = NaN;
        intf_meta(m).d_REFL_d2    = NaN;
        intf_meta(m).d_REFL_total = NaN;
    end

    othersIDs_V(m) = I(m).SensorID;
end

%% ── Section C: Waveform Cache Initialisation ─────────────────────────────
%   The waveform cache is initialised on the first call.  It persists
%   across CPI micro-steps via ctx and avoids regenerating LO signals
%   for static schemes (where F0 does not change between CPIs).

if ~isfield(ctx, 'waveformCache') || isempty(ctx.waveformCache)
    ctx.waveformCache = struct();
end

%% ── Section D: Victim LO Signal Generation ───────────────────────────────
%   Generate the victim's time–frequency LO signal.
%   For static schemes (ARIFOMA_scheme_is_dynamic = false), the signal is
%   cached and reused in subsequent CPIs.  Dynamic schemes always regenerate.

keyV      = sprintf('S%d', V.SensorID);
useCacheV = ~ARIFOMA_scheme_is_dynamic(V.Scheme);

if useCacheV && isfield(ctx.waveformCache, keyV)
    % Retrieve cached signal — no regeneration needed.
    base = ctx.waveformCache.(keyV);
else
    % Generate a fresh LO signal and cache if scheme is static.
    [t_base, f_base, idx_base] = ARIFOMA_generate_LO_Signal_by_Scheme( ...
        V.Scheme, V.Waveform, fv_base, Bv, Tv, Nch_V, ...
        T_flyback_default_V, T_flyback_stepped_V, T_returnfly_V, ...
        chirpsPerStep_V, freqIncr_V, dir_V, +1, ...
        bandMin_V, bandMax_V, ...
        Nseg_V, PermPattern_V);

    base = struct('t', t_base, 'f', f_base, 'idx', idx_base);

    if useCacheV
        ctx.waveformCache.(keyV) = base;
    end
end

% Apply a random start-time offset to model unsynchronised transmission.
% Applied twice to increase randomisation of the victim's timing.
t_v   = base.t + rand * cfg.randstartTime;
t_v   = t_v    + rand * cfg.randstartTime;
f_v   = base.f;
idx_v = base.idx;
Nch_V = max(idx_v);   % update Nch_V from the generated signal

%% ── Section E: Interferer LO Signal Generation ───────────────────────────
%   Generate LO signals for all interferers in this corner view.
%   Caching logic mirrors Section D.

interferer_t_all   = cell(1, nI);
interferer_f_all   = cell(1, nI);
interferer_idx_all = cell(1, nI);

for j = 1:nI

    keyI      = sprintf('S%d', I(j).SensorID);
    schemeI   = I(j).Scheme;
    useCacheI = ~ARIFOMA_scheme_is_dynamic(schemeI);

    if useCacheI && isfield(ctx.waveformCache, keyI)
        base = ctx.waveformCache.(keyI);
    else
        [t_base, f_base, idx_base] = ARIFOMA_generate_LO_Signal_by_Scheme( ...
            schemeI, I(j).Waveform, fi_base(j), ...
            I(j).B, I(j).T, I(j).Nch, ...
            I(j).T_flyback_default, I(j).T_flyback_stepped, I(j).T_returnfly, ...
            I(j).ChirpsPerStep, I(j).FreqIncrHz, I(j).Dir, +1, ...
            I(j).bandMin, I(j).bandMax, ...
            I(j).Nseg, I(j).PermPattern);

        base = struct('t', t_base, 'f', f_base, 'idx', idx_base);

        if useCacheI
            ctx.waveformCache.(keyI) = base;
        end
    end

    % Apply random start-time offset per interferer.
    interferer_t_all{j}   = base.t + rand * cfg.randstartTime;
    interferer_f_all{j}   = base.f;
    interferer_idx_all{j} = base.idx;

end   % interferer generation loop

%% ── Section F: Analytic Interference Detection and FOM Accumulation ──────
%   For each interferer, run the analytic detector against the victim
%   waveform and accumulate per-chirp and per-path-type metrics into bins.
%
%   collMask(k,j) = true if victim chirp k overlaps with interferer j.
%   The per-CPI interference ratio counts chirps where any column is true.

chIDs = unique(idx_v(:).', 'stable');
Nch   = numel(chIDs);
M     = nI;

% Per-chirp interferer count arrays (one element per victim chirp).
k_total = zeros(1, Nch);
k_LOS   = zeros(1, Nch);
k_REFL  = zeros(1, Nch);

% Collision mask: rows = victim chirps, columns = interferers.
collMask = false(Nch, M);

% Accumulation bins for interference durations and distances.
bins.LOS.dur   = [];
bins.LOS.dist  = [];
bins.REFL.d1   = [];
bins.REFL.d2   = [];
bins.REFL.dur  = [];
bins.REFL.dist = [];

% Global interference interval start/end times for spectrogram overlay.
t0_all = [];
t1_all = [];

for j = 1:M

    tj = interferer_t_all{j};
    fj = interferer_f_all{j};
    ij = interferer_idx_all{j};

    if isempty(tj),  continue;  end

    % ── Analytic detector: solve |f_V(t) - f_I(t)| <= BADC ───────────────
    [dt_j, ~, ts_j, te_j, mask_j, ~, ~, ~] = ...
        ARIFOMA_analytic_detector( ...
            t_v, f_v, idx_v, tj, fj, ij, ...
            cfg.BADC_Hz, cfg.boolIgnoreFlybackInterferer);

    if isempty(mask_j),  continue;  end

    mask_j            = logical(mask_j(:).');
    collMask(:, j)    = mask_j;
    k_total           = k_total + mask_j;

    % ── Bin by path type ──────────────────────────────────────────────────
    if isfinite(intf_meta(j).d_REFL_total) && intf_meta(j).d_REFL_total > 0
        % Reflected path contribution.
        k_REFL          = k_REFL + mask_j;
        bins.REFL.dur   = [bins.REFL.dur,  dt_j(:).'];
        bins.REFL.dist  = [bins.REFL.dist, repmat(intf_meta(j).d_REFL_total, 1, numel(dt_j))];
        bins.REFL.d1    = [bins.REFL.d1,   repmat(intf_meta(j).d_REFL_d1,    1, numel(dt_j))];
        bins.REFL.d2    = [bins.REFL.d2,   repmat(intf_meta(j).d_REFL_d2,    1, numel(dt_j))];
    elseif isfinite(intf_meta(j).d_LOS) && intf_meta(j).d_LOS > 0
        % LOS path contribution.
        k_LOS           = k_LOS + mask_j;
        bins.LOS.dur    = [bins.LOS.dur,  dt_j(:).'];
        bins.LOS.dist   = [bins.LOS.dist, repmat(intf_meta(j).d_LOS, 1, numel(dt_j))];
    end

    % Accumulate global interval times for spectrogram overlay.
    if ~isempty(ts_j)
        t0_all = [t0_all; ts_j(:)];   %#ok<AGROW>
        t1_all = [t1_all; te_j(:)];   %#ok<AGROW>
    end

end   % interferer detection loop

% ── CPI-level interference ratio ─────────────────────────────────────────
%   Count chirps where at least one interferer causes overlap.
total_interfered_chirps = sum(any(collMask, 2));
P_coll_per_n            = total_interfered_chirps / Nch;

%% ── Section G: Pack Output Struct ───────────────────────────────────────
%   All FOMs and the updated context are packed into a flat output struct.
%   Interference durations are normalised by the victim chirp duration Tv.

outS = struct();
outS.ctx = ctx;

% CPI-level interference ratio FOM.
outS.P_coll_per_n           = P_coll_per_n;
outS.total_interfered_chirps = total_interfered_chirps;

% Normalised interference durations (T_I / T_chirp).
outS.all_interference_durations      = [bins.LOS.dur, bins.REFL.dur] / Tv;
outS.all_interference_durations_LOS  = bins.LOS.dur  / Tv;
outS.all_interference_durations_REFL = bins.REFL.dur / Tv;

% Path distances per interference occurrence.
outS.all_interference_dist_LOS       = bins.LOS.dist;
outS.all_interference_dist_REFL      = bins.REFL.dist;
outS.all_interference_dist_REFL_d1   = bins.REFL.d1;
outS.all_interference_dist_REFL_d2   = bins.REFL.d2;

% Per-chirp interferer counts.
outS.all_nIntf_per_chirp_total = k_total;
outS.all_nIntf_per_chirp_LOS   = k_LOS;
outS.all_nIntf_per_chirp_REFL  = k_REFL;

% BATS centroid estimate — not computed in the conventional strategy.
outS.all_r_hat_sub = NaN;

%% ── Section H: Time–Frequency Spectrogram Visualisation ─────────────────
%   Renders the victim and interferer LO signals as a time–frequency
%   spectrogram overlay on the corner tile axes.
%   Rendered only when cfg.plotScenarioEnd = true.
%
%   Colour convention:
%     Blue  — victim chirp
%     Red   — LOS interferer (semi-transparent)
%     Grey  — reflected interferer (semi-transparent)
%     Blue shaded patch — interference region

if cfg.plotScenarioEnd

    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'on');

    % ── Victim waveform ───────────────────────────────────────────────────
    plot(ax, t_v*1e3, f_v*1e-9, 'b', ...
        'DisplayName','Victim $\nu$', 'LineWidth',2);

    % ── Interferer waveforms (only when interferers exist) ────────────────
    if ~strcmpi(schemeInterferer, "NOTHING")

        alpha_ifr = 0.30;                     % transparency for all interferers
        col_LOS   = [0.90, 0.30, 0.30];       % red  — LOS interferer
        col_REFL  = [0.65, 0.65, 0.65];       % grey — reflected interferer

        for j = 1:nI
            % Select colour by path type.
            isRefl = isfield(I(j), 'min_d2_m') && ...
                     ~isempty(I(j).min_d2_m)   && ...
                     isfinite(I(j).min_d2_m)   && ...
                     I(j).min_d2_m > 0;

            col   = col_LOS  * ~isRefl + col_REFL * isRefl;
            label = sprintf('Interferer $\\gamma$ %d', j);

            plot(ax, interferer_t_all{j}*1e3, interferer_f_all{j}*1e-9, ...
                'Color',    [col, alpha_ifr], ...   % RGBA: HG2 supported
                'LineStyle', '-', ...
                'LineWidth', 2.0, ...
                'DisplayName', label);
        end

        % ── Interference region shading ───────────────────────────────────
        %   Blue patches mark all detected interference intervals.
        %   Only the first patch gets a legend entry to avoid duplication.
        addedInterfLegend = false;

        for k = 1:numel(t0_all)
            x_patch = [t0_all(k), t1_all(k), t1_all(k), t0_all(k)] * 1e3;
            y_patch = [bandMin_V, bandMin_V, bandMax_V, bandMax_V]   * 1e-9;

            if ~addedInterfLegend
                fill(ax, x_patch, y_patch, 'b', ...
                    'FaceAlpha',       0.4, ...
                    'EdgeColor',       'none', ...
                    'DisplayName',     'Interference');
                addedInterfLegend = true;
            else
                fill(ax, x_patch, y_patch, 'b', ...
                    'FaceAlpha',        0.4, ...
                    'EdgeColor',        'none', ...
                    'HandleVisibility', 'off', ...
                    'Tag',              'interfRegion');
            end
        end

    end   % interferer rendering

    % ── Axes formatting ───────────────────────────────────────────────────
    ylim(ax, [bandMin_V*1e-9, bandMax_V*1e-9]);
    xlim(ax, [0, cfg.CPI_duration_s*1e3]);
    xlabel(ax, 'Time (ms)');
    ylabel(ax, 'Frequency (GHz)');
    title(ax, sprintf('Corner Radar %s | $n_\\mathrm{F}=$%d', cIndex, numCPI));
    legend(ax, 'show', 'Location','eastoutside');
    grid(ax, 'on');

end   % plotScenarioEnd

end   % ARIFOMA_conventional_scene_runCPI

% =========================================================================
%% ── Local Helper Functions ───────────────────────────────────────────────
% =========================================================================

% -------------------------------------------------------------------------
function tf = ARIFOMA_scheme_is_dynamic(scheme)
% ARIFOMA_SCHEME_IS_DYNAMIC
%   Return true if a mitigation scheme requires fresh LO signal generation
%   at every CPI micro-step (i.e. the carrier frequency may change).
%
%   Dynamic schemes hop F0 per CPI or per PRI and must not use the
%   waveform cache.  Static schemes keep F0 fixed across CPIs and can
%   safely reuse a cached LO signal.
%
%   Dynamic schemes (no cache):
%     FAR                           — frequency-agile radar (per-PRI hop)
%     RANDOM_FH_PER_CPI             — random hop at each CPI
%     RANDOM_FH_PER_CPI_INTERFERED  — conditional random hop at each CPI
%     BATSINSPIRED_FH_PER_CPI       — BATS-directed hop at each CPI
%
%   INPUT
%     scheme  (string) — scheme identifier string
%
%   OUTPUT
%     tf  (logical) — true if the scheme is dynamic (no caching allowed)

    tf = ismember(scheme, [ ...
        "FAR", ...
        "RANDOM_FH_PER_CPI", ...
        "RANDOM_FH_PER_CPI_INTERFERED", ...
        "BATSINSPIRED_FH_PER_CPI"]);
end