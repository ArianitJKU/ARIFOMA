function outS = ARIFOMA_batsinspired_scene_runCPI( ...
        ax, numCPI, cfg, cIndex, ctx, V, I, ...
        fv_base, fi_base, schemeInterferer)
% ARIFOMA_BATSINSPIRED_SCENE_RUNCPI
%   Execute one CPI micro-step for a single victim corner using the
%   BATS-inspired frequency-hopping interference analysis pipeline.
%
% PURPOSE
% -------
%   This function extends ARIFOMA_conventional_scene_runCPI with the
%   additional BATS (bat-inspired) observables required by the
%   BATSINSPIRED_FH_PER_CPI mitigation scheme.  The key difference from
%   the conventional pipeline is the use of ARIFOMA_analytic_detector_batsinspired
%   instead of ARIFOMA_analytic_detector, which additionally returns:
%
%     r_hat_j     — per-interferer interference centroid estimate [0, 1]
%     f_int_min_j — minimum victim RF frequency that was interfered [Hz]
%     f_int_max_j — maximum victim RF frequency that was interfered [Hz]
%
%   These observables are accumulated across all interferers using a
%   duration-weighted mean and stored in outS for use by
%   ARIFOMA_apply_interferer_schemes in the next CPI micro-step.
%
% BATS OBSERVABLE ACCUMULATION
% ----------------------------
%   For each interferer j that produces at least one interference interval:
%
%     r_num_all += r_hat_j * sum(dt_j)    (duration-weighted numerator)
%     r_den_all += sum(dt_j)              (total interfered duration)
%
%   The CPI-level centroid estimate is then:
%
%     outS.all_r_hat_sub = r_num_all / r_den_all
%
%   This is the continuous-time analogue of the discrete BATS estimator
%   from [Bechter et al., IEEE MTT-S ICMIM 2016].
%
%   f_int_min and f_int_max track the envelope of the interfered victim
%   RF frequencies across all interferers, enabling directional band shift.
%
% DIFFERENCE FROM CONVENTIONAL PIPELINE
% --------------------------------------
%   Section F (interference detection loop) is the only section that
%   differs from ARIFOMA_conventional_scene_runCPI:
%     - Uses ARIFOMA_analytic_detector_batsinspired instead of
%       ARIFOMA_analytic_detector
%     - Accumulates r_hat, f_int_min, f_int_max per interferer
%   All other sections (waveform generation, caching, FOM binning,
%   visualisation) are identical to the conventional pipeline.
%
% INPUTS
%   ax                (axes)   — target corner tile axes handle
%   numCPI            (int)    — current CPI micro-step index (1-based)
%   cfg               (struct) — global configuration struct; required fields:
%                                  cfg.randstartTime           [s]
%                                  cfg.BADC_Hz                 [Hz]
%                                  cfg.boolIgnoreFlybackInterferer (logical)
%                                  cfg.CPI_duration_s          [s]
%                                  cfg.plotScenarioEnd         (logical)
%   cIndex            (int or string) — corner identifier for tile title
%   ctx               (struct) — current corner context; modified fields:
%                                  ctx.waveformCache
%   V                 (struct) — victim sensor parameter row
%   I                 (struct array, 1×nI) — interferer parameter rows;
%                                  each must contain min_d1_m and min_d2_m
%   fv_base           (double) — victim carrier frequency F0 [Hz]
%   fi_base           (double, nI×1) — interferer carrier frequencies [Hz]
%   schemeInterferer  (string) — active interferer scheme identifier
%
% OUTPUT
%   outS  (struct) — per-CPI result struct with all fields from
%                    ARIFOMA_conventional_scene_runCPI, plus:
%
%     BATS observables
%       outS.all_r_hat_sub   — duration-weighted CPI-level centroid estimate
%                              r_hat ∈ [0,1]; NaN if no interference detected
%       outS.f_int_min       — minimum interfered victim RF frequency [Hz];
%                              NaN if no interference detected
%       outS.f_int_max       — maximum interfered victim RF frequency [Hz];
%                              NaN if no interference detected
%
% REFERENCES
%   [1] J. Bechter, C. Sippel, and C. Waldschmidt, "Bats-inspired frequency
%       hopping for mitigation of interference between automotive radars,"
%       IEEE MTT-S ICMIM, 2016.
%
% SEE ALSO
%   ARIFOMA_conventional_scene_runCPI,
%   ARIFOMA_analytic_detector_batsinspired,
%   ARIFOMA_apply_interferer_schemes,
%   ARIFOMA_generate_LO_Signal_by_Scheme

% =========================================================================

%% ── Section A: Victim Waveform Parameters ────────────────────────────────

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
%   Classify each interferer's dominant path as LOS or reflected and
%   store the relevant distances for FOM binning in Section F.

nI        = numel(I);
intf_meta = repmat(struct( ...
    'd_LOS',NaN, 'd_REFL_total',NaN, ...
    'd_REFL_d1',NaN, 'd_REFL_d2',NaN), 1, nI);

othersIDs_V = zeros(1, nI);

for m = 1:nI
    d1 = double(I(m).min_d1_m);
    d2 = double(I(m).min_d2_m);

    if isfinite(d2) && d2 > 0
        % Reflected path: both legs defined.
        intf_meta(m).d_LOS        = NaN;
        intf_meta(m).d_REFL_d1    = d1;
        intf_meta(m).d_REFL_d2    = d2;
        intf_meta(m).d_REFL_total = d1 + d2;
    else
        % LOS path: only first leg defined.
        intf_meta(m).d_LOS        = d1;
        intf_meta(m).d_REFL_d1    = NaN;
        intf_meta(m).d_REFL_d2    = NaN;
        intf_meta(m).d_REFL_total = NaN;
    end

    othersIDs_V(m) = I(m).SensorID;
end

%% ── Section C: Waveform Cache Initialisation ─────────────────────────────
if ~isfield(ctx, 'waveformCache') || isempty(ctx.waveformCache)
    ctx.waveformCache = struct();
end

%% ── Section D: Victim LO Signal Generation (Cached) ─────────────────────
%   Static schemes reuse the cached LO signal; dynamic schemes regenerate.

keyV      = sprintf('S%d', V.SensorID);
useCacheV = ~ARIFOMA_scheme_is_dynamic(V.Scheme);

if useCacheV && isfield(ctx.waveformCache, keyV)
    base = ctx.waveformCache.(keyV);
else
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

% Apply random start-time offset (twice) to model unsynchronised TX.
t_v   = base.t + rand * cfg.randstartTime;
t_v   = t_v    + rand * cfg.randstartTime;
f_v   = base.f;
idx_v = base.idx;
Nch_V = max(idx_v);

%% ── Section E: Interferer LO Signal Generation (Cached) ─────────────────
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

    interferer_t_all{j}   = base.t + rand * cfg.randstartTime;
    interferer_f_all{j}   = base.f;
    interferer_idx_all{j} = base.idx;
end

%% ── Section F: BATS-Extended Interference Detection and FOM Accumulation ─
%   This section differs from the conventional pipeline by calling
%   ARIFOMA_analytic_detector_batsinspired, which additionally returns the
%   per-interferer interference centroid r_hat_j and the interfered RF
%   frequency range [f_int_min_j, f_int_max_j].
%
%   These are accumulated across all interferers using duration weighting
%   to produce the CPI-level BATS observables stored in outS.

chIDs = unique(idx_v(:).', 'stable');
Nch   = numel(chIDs);
M     = nI;

k_total  = zeros(1, Nch);
k_LOS    = zeros(1, Nch);
k_REFL   = zeros(1, Nch);
collMask = false(Nch, M);

bins.LOS.dur   = [];
bins.LOS.dist  = [];
bins.REFL.d1   = [];
bins.REFL.d2   = [];
bins.REFL.dur  = [];
bins.REFL.dist = [];

t0_all = [];
t1_all = [];

% ── BATS accumulator initialisation ──────────────────────────────────────
%   r_num_all / r_den_all accumulate the duration-weighted centroid sum.
%   f_int_min/max_all track the envelope of interfered victim RF frequencies.
r_num_all     = 0;
r_den_all     = 0;
f_int_min_all = NaN;
f_int_max_all = NaN;

for j = 1:M

    tj = interferer_t_all{j};
    fj = interferer_f_all{j};
    ij = interferer_idx_all{j};
    if isempty(tj),  continue;  end

    % ── BATS-extended analytic detector ───────────────────────────────────
    %   Returns the standard FOMs plus r_hat_j, f_int_min_j, f_int_max_j.
    [dt_j, ~, ts_j, te_j, mask_j, ~, ~, ~, r_hat_j, f_int_min_j, f_int_max_j] = ...
        ARIFOMA_analytic_detector_batsinspired( ...
            t_v, f_v, idx_v, tj, fj, ij, ...
            cfg.BADC_Hz, cfg.boolIgnoreFlybackInterferer);

    % ── Accumulate BATS centroid and frequency range ───────────────────────
    if isfinite(r_hat_j) && ~isempty(dt_j)
        dur_j     = sum(dt_j);
        r_num_all = r_num_all + r_hat_j * dur_j;
        r_den_all = r_den_all + dur_j;

        % Expand the interfered frequency envelope.
        if isfinite(f_int_min_j)
            if isnan(f_int_min_all)
                f_int_min_all = f_int_min_j;
                f_int_max_all = f_int_max_j;
            else
                f_int_min_all = min(f_int_min_all, f_int_min_j);
                f_int_max_all = max(f_int_max_all, f_int_max_j);
            end
        end
    end

    if isempty(mask_j),  continue;  end

    % ── Standard FOM accumulation (identical to conventional pipeline) ────
    mask_j         = logical(mask_j(:).');
    collMask(:, j) = mask_j;
    k_total        = k_total + mask_j;

    if isfinite(intf_meta(j).d_REFL_total) && intf_meta(j).d_REFL_total > 0
        k_REFL          = k_REFL + mask_j;
        bins.REFL.dur   = [bins.REFL.dur,  dt_j(:).'];
        bins.REFL.dist  = [bins.REFL.dist, repmat(intf_meta(j).d_REFL_total, 1, numel(dt_j))];
        bins.REFL.d1    = [bins.REFL.d1,   repmat(intf_meta(j).d_REFL_d1,    1, numel(dt_j))];
        bins.REFL.d2    = [bins.REFL.d2,   repmat(intf_meta(j).d_REFL_d2,    1, numel(dt_j))];
    elseif isfinite(intf_meta(j).d_LOS) && intf_meta(j).d_LOS > 0
        k_LOS           = k_LOS + mask_j;
        bins.LOS.dur    = [bins.LOS.dur,  dt_j(:).'];
        bins.LOS.dist   = [bins.LOS.dist, repmat(intf_meta(j).d_LOS, 1, numel(dt_j))];
    end

    if ~isempty(ts_j)
        t0_all = [t0_all; ts_j(:)];   %#ok<AGROW>
        t1_all = [t1_all; te_j(:)];   %#ok<AGROW>
    end

end   % interferer detection loop

% ── CPI-level interference ratio ─────────────────────────────────────────
total_interfered_chirps = sum(any(collMask, 2));
P_coll_per_n            = total_interfered_chirps / Nch;

%% ── Section G: Pack Output Struct ───────────────────────────────────────
outS = struct();
outS.ctx = ctx;

% CPI-level interference ratio FOM.
outS.P_coll_per_n            = P_coll_per_n;
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

% ── BATS observables ──────────────────────────────────────────────────────
%   Duration-weighted CPI-level centroid estimate.
%   NaN if no interference was detected across all interferers.
if r_den_all > 0
    outS.all_r_hat_sub = r_num_all / r_den_all;
else
    outS.all_r_hat_sub = NaN;
end

outS.f_int_min = f_int_min_all;
outS.f_int_max = f_int_max_all;

%% ── Section H: Time–Frequency Spectrogram Visualisation ─────────────────
%   Identical to ARIFOMA_conventional_scene_runCPI.
%   Rendered only when cfg.plotScenarioEnd = true.

if cfg.plotScenarioEnd

    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'on');

    % Victim waveform.
    plot(ax, t_v*1e3, f_v*1e-9, 'b', ...
        'DisplayName','Victim $\nu$', 'LineWidth',2);

    if ~strcmpi(schemeInterferer, "NOTHING")

        alpha_ifr = 0.30;
        col_LOS   = [0.90, 0.30, 0.30];   % red  — LOS
        col_REFL  = [0.65, 0.65, 0.65];   % grey — reflected

        for j = 1:nI
            isRefl = isfield(I(j), 'min_d2_m') && ...
                     ~isempty(I(j).min_d2_m)   && ...
                     isfinite(I(j).min_d2_m)   && ...
                     I(j).min_d2_m > 0;

            col   = col_LOS * ~isRefl + col_REFL * isRefl;
            label = sprintf('Interferer $\\gamma$ %d', j);

            plot(ax, interferer_t_all{j}*1e3, interferer_f_all{j}*1e-9, ...
                'Color',     [col, alpha_ifr], ...
                'LineStyle', '-', ...
                'LineWidth', 2.0, ...
                'DisplayName', label);
        end

        % Interference region patches — single legend entry.
        addedInterfLegend = false;
        for k = 1:numel(t0_all)
            x_patch = [t0_all(k), t1_all(k), t1_all(k), t0_all(k)] * 1e3;
            y_patch = [bandMin_V, bandMin_V, bandMax_V, bandMax_V]   * 1e-9;

            if ~addedInterfLegend
                fill(ax, x_patch, y_patch, 'b', ...
                    'FaceAlpha',   0.4, ...
                    'EdgeColor',   'none', ...
                    'DisplayName', 'Interference');
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

    ylim(ax, [bandMin_V*1e-9, bandMax_V*1e-9]);
    xlim(ax, [0, cfg.CPI_duration_s*1e3]);
    xlabel(ax, 'Time (ms)');
    ylabel(ax, 'Frequency (GHz)');
    title(ax, sprintf('Corner Radar %s | $n_\\mathrm{F}=$%d', cIndex, numCPI));
    legend(ax, 'show', 'Location','eastoutside');
    grid(ax, 'on');

end   % plotScenarioEnd

end   % ARIFOMA_batsinspired_scene_runCPI

% =========================================================================
%% ── Local Helper Functions ───────────────────────────────────────────────
% =========================================================================

% -------------------------------------------------------------------------
function tf = ARIFOMA_scheme_is_dynamic(scheme)
% ARIFOMA_SCHEME_IS_DYNAMIC
%   Return true if a mitigation scheme requires fresh LO signal generation
%   at every CPI micro-step.
%
%   Dynamic schemes (no cache):
%     FAR, RANDOM_FH_PER_CPI, RANDOM_FH_PER_CPI_INTERFERED,
%     BATSINSPIRED_FH_PER_CPI
%
%   INPUT
%     scheme  (string) — scheme identifier string
%
%   OUTPUT
%     tf  (logical) — true if the scheme is dynamic

    tf = ismember(scheme, [ ...
        "FAR", ...
        "RANDOM_FH_PER_CPI", ...
        "RANDOM_FH_PER_CPI_INTERFERED", ...
        "BATSINSPIRED_FH_PER_CPI"]);
end