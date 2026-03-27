function engine = ARIFOMA_build_engine_inputs( ...
        ax, corner, pathByCorner, cfg, ...
        numCPI, numberofCPIsperAcquisition, tIndex, t0_cur)
% ARIFOMA_BUILD_ENGINE_INPUTS
%   Assemble the engine input bundle consumed by ARIFOMA_run_strategy.
%
% PURPOSE
% -------
%   This function is the final preparation stage before strategy dispatch.
%   It collects all objects needed for waveform generation and interference
%   analysis into a single flat struct (engine) so that ARIFOMA_run_strategy
%   and its sub-functions receive a consistent, self-contained bundle
%   regardless of which corner or CPI micro-step is being processed.
%
%   Responsibilities:
%     1. Extract per-corner LOS and reflected path distances from
%        pathByCorner.
%     2. Normalise the corner state struct to guarantee consistent field
%        layout for victim and all interferers.
%     3. Initialise fresh per-CPI accumulators and log structs.
%     4. Derive the chirp index array iv_base from the absolute time grid
%        if it was not already computed upstream.
%     5. Pack all objects into engine.
%
% PIPELINE POSITION
% -----------------
%   Called from runCPI_perCorner (main_sim.m) as Stage 4:
%
%     Stage 3                         Stage 4                    Stage 5
%   ARIFOMA_apply_all_schemes  →  ARIFOMA_build_engine_inputs  →  ARIFOMA_run_strategy
%
% INPUTS
%   ax                         (axes)   — target tile axes handle for
%                                         in-loop visualisation
%   corner                     (struct) — updated corner bundle from
%                                         ARIFOMA_apply_all_schemes;
%                                         must contain: V, I, ctx, cn,
%                                         fv_base (optional), iv_base (optional)
%   pathByCorner               (struct) — per-corner path distance struct
%                                         indexed by corner name (e.g. FL);
%                                         each entry may contain:
%                                           dist_LOS_m, dist_REFL_total_m,
%                                           dist_REFL_d1_m, dist_REFL_d2_m
%   cfg                        (struct) — global configuration struct
%   numCPI                     (int)    — current CPI micro-step index (1-based)
%   numberofCPIsperAcquisition (int)    — total CPI micro-steps per snapshot
%   tIndex                     (int)    — current acquisition time index (1-based)
%   t0_cur                     (double) — current acquisition timestamp [s]
%
% OUTPUT
%   engine  (struct) — self-contained engine input bundle with fields:
%
%     Visualisation
%       engine.ax                          axes handle
%
%     Configuration and timing
%       engine.cfg                         global config struct
%       engine.t0_cur                      acquisition timestamp [s]
%       engine.numCPI                      CPI micro-step index
%       engine.numberofCPIsperAcquisition  total micro-steps per snapshot
%       engine.tIndex                      acquisition time index
%       engine.t_base                      absolute time grid [s]
%
%     Scene objects
%       engine.corner                      full corner bundle
%       engine.ctx                         corner context struct
%       engine.state                       normalised state struct
%
%     Waveform indices
%       engine.fv_base                     victim carrier frequency [Hz]
%       engine.iv_base                     chirp index array (NaN during flyback)
%
%     Path distances
%       engine.dist.LOS                    LOS path distances [m]
%       engine.dist.REFL_total             total reflected path distances [m]
%       engine.dist.REFL_d1                first reflected leg distances [m]
%       engine.dist.REFL_d2                second reflected leg distances [m]
%
%     Accumulators
%       engine.acc                         fresh per-CPI accumulator struct
%
% NOTES
% -----
%   - engine.t_base is constructed as t0_cur + tIndex(:).  When tIndex is
%     a scalar this yields a single-element time grid; when tIndex is a
%     vector it yields the full absolute time axis for the snapshot.
%   - iv_base is derived from make_idx_from_time if corner.iv_base is absent
%     or empty.  The derivation uses V.T (chirp duration) and V.T_fly
%     (flyback duration, defaulting to 0 if absent).
%   - Flyback samples are marked NaN in iv_base so that chirp masks in
%     downstream detectors correctly exclude non-TX intervals.
%
% SEE ALSO
%   ARIFOMA_apply_all_schemes,
%   ARIFOMA_run_strategy,
%   ARIFOMA_init_accumulators,
%   main_sim

% =========================================================================

%% ── Unpack Corner Bundle ─────────────────────────────────────────────────
ctx = corner.ctx;
V   = corner.V;
I   = corner.I;   %#ok<NASGU>  retained for clarity; accessed via corner

%% ── Step 1: Extract Per-Corner Path Distances ────────────────────────────
%   Pull LOS and reflected path distance arrays from pathByCorner for this
%   corner.  Missing fields are returned as empty arrays (safe fallback).

[dist_LOS_m, dist_REFL_total_m, dist_REFL_d1_m, dist_REFL_d2_m] = ...
    extract_corner_distances(pathByCorner, corner.cn);

%% ── Step 2: Normalise State Struct ──────────────────────────────────────
%   Ensure state.victim and state.I(1:nI) exist with consistent layouts.
%   Also propagate the victim sweep bandwidth and chirp duration into state
%   so downstream functions do not need to re-read V.

state          = normalize_state(ctx.state, numel(I));
state.victim.B = V.B;   % sweep bandwidth [Hz]
state.victim.T = V.T;   % chirp duration [s]

%% ── Step 4: Pack Engine Bundle ───────────────────────────────────────────
%   All objects are assembled into a single flat struct.  Flat packing
%   avoids deep nesting that would require callers to know the internal
%   structure of corner or ctx.

engine.ax                          = ax;
engine.cfg                         = cfg;
engine.corner                      = corner;
engine.ctx                         = ctx;
engine.state                       = state;
engine.t0_cur                      = t0_cur;
engine.numCPI                      = numCPI;
engine.numberofCPIsperAcquisition  = numberofCPIsperAcquisition;
engine.tIndex                      = tIndex;

% Absolute time grid for the current acquisition snapshot.
% Strategies should use this grid rather than constructing their own.
engine.t_base = t0_cur + tIndex(:);

%% ── Step 5: Victim Carrier Frequency ────────────────────────────────────
%   Copy fv_base from the corner bundle if already computed by
%   ARIFOMA_apply_all_schemes.

if isfield(corner, 'fv_base')
    engine.fv_base = corner.fv_base;
end

%% ── Step 6: Chirp Index Array (iv_base) ─────────────────────────────────
%   iv_base maps each time sample in t_base to its chirp index (1, 2, …).
%   Flyback samples are marked NaN so detector chirp masks exclude them.
%
%   If iv_base was passed in via corner (e.g. pre-computed upstream), use
%   it directly.  Otherwise derive it from the victim's chirp geometry.

if isfield(corner, 'iv_base') && ~isempty(corner.iv_base)
    engine.iv_base = corner.iv_base;
else
    % Derive iv_base from V.T (chirp duration) and V.T_fly (flyback).
    T    = V.T;
    Tfly = 0;
    if isfield(V, 'T_fly'),  Tfly = V.T_fly;  end

    engine.iv_base = make_idx_from_time(engine.t_base, t0_cur, T, Tfly);
end

%% ── Step 7: Path Distance Sub-Struct ────────────────────────────────────
engine.dist = struct( ...
    'LOS',        dist_LOS_m, ...
    'REFL_total', dist_REFL_total_m, ...
    'REFL_d1',    dist_REFL_d1_m, ...
    'REFL_d2',    dist_REFL_d2_m);

end   % ARIFOMA_build_engine_inputs

% =========================================================================
%% ── Local Helper Functions ───────────────────────────────────────────────
% =========================================================================

% -------------------------------------------------------------------------
function [dist_LOS_m, dist_REFL_total_m, dist_REFL_d1_m, dist_REFL_d2_m] = ...
        extract_corner_distances(pathByCorner, cn)
% EXTRACT_CORNER_DISTANCES
%   Extract per-corner path distance arrays from the pathByCorner struct.
%
% PURPOSE
% -------
%   Retrieves the LOS and reflected path distance arrays for corner cn
%   from pathByCorner (produced by ARIFOMA_prepare_scene_views).  If the
%   corner field or any distance sub-field is absent, the corresponding
%   output is returned as an empty array.  This makes the function safe to
%   call even when path data is partially available.
%
% INPUTS
%   pathByCorner  (struct) — per-corner path distance struct indexed by
%                            corner name; each entry may contain:
%                              dist_LOS_m, dist_REFL_total_m,
%                              dist_REFL_d1_m, dist_REFL_d2_m
%   cn            (char)  — corner name string (e.g. 'FL', 'RR')
%
% OUTPUTS
%   dist_LOS_m        (double array) — LOS path distances [m]; [] if absent
%   dist_REFL_total_m (double array) — total reflected path distances [m]
%   dist_REFL_d1_m    (double array) — first reflected leg distances [m]
%   dist_REFL_d2_m    (double array) — second reflected leg distances [m]

    dist_LOS_m        = [];
    dist_REFL_total_m = [];
    dist_REFL_d1_m    = [];
    dist_REFL_d2_m    = [];

    try
        P = pathByCorner.(cn);

        if isfield(P, 'dist_LOS_m'),        dist_LOS_m        = P.dist_LOS_m;        end
        if isfield(P, 'dist_REFL_total_m'), dist_REFL_total_m = P.dist_REFL_total_m; end
        if isfield(P, 'dist_REFL_d1_m'),    dist_REFL_d1_m    = P.dist_REFL_d1_m;    end
        if isfield(P, 'dist_REFL_d2_m'),    dist_REFL_d2_m    = P.dist_REFL_d2_m;    end
    catch
        % Corner field absent or pathByCorner is empty — leave outputs as [].
    end
end

% -------------------------------------------------------------------------
function state = normalize_state(state, nI)
% NORMALIZE_STATE
%   Ensure the state struct has consistent victim and interferer sub-structs.
%
% PURPOSE
% -------
%   The state struct is built incrementally across CPI micro-steps and may
%   not always have fully populated sub-fields on the first call or after
%   the number of interferers changes between snapshots.  This function
%   guarantees that state.victim and state.I(1:nI) always exist with at
%   least a minimal valid layout before engine packing.
%
% INPUTS
%   state  (struct) — current state struct from ctx.state
%   nI     (int)    — expected number of interferer entries in state.I
%
% OUTPUT
%   state  (struct) — normalised state with guaranteed fields:
%                       state.victim  — at minimum {has: false}
%                       state.I       — struct array of length nI,
%                                       each entry at minimum {has: false}

    % Ensure victim sub-struct exists.
    if ~isfield(state, 'victim') || isempty(state.victim)
        state.victim = struct('has', false);
    end

    % Ensure interferer array exists with correct length.
    if ~isfield(state, 'I') || isempty(state.I)
        state.I = repmat(struct('has', false), 1, nI);
    elseif numel(state.I) < nI
        % Extend array if new interferers appeared since last snapshot.
        state.I(end+1:nI) = repmat(struct('has', false), 1, nI - numel(state.I));
    elseif numel(state.I) > nI
        % Trim array if interferers dropped out since last snapshot.
        state.I = state.I(1:nI);
    end
end

% -------------------------------------------------------------------------
function idx = make_idx_from_time(t_abs, t0, T, Tfly)
% MAKE_IDX_FROM_TIME
%   Build a chirp index array from an absolute time grid.
%
% PURPOSE
% -------
%   Maps each time sample in t_abs to its 1-based chirp index within the
%   CPI waveform.  Samples that fall during the flyback interval (i.e.
%   after the active ramp but before the next chirp starts) are assigned
%   NaN so that downstream chirp masks correctly exclude non-TX intervals.
%
% ALGORITHM
% ---------
%   1. Compute tau = t_abs - t0  (time relative to CPI start).
%   2. Compute the chirp period:  Tper = T + Tfly.
%   3. Assign chirp index:        k = floor(tau / Tper) + 1  (1-based).
%   4. Compute phase within PRI:  ph = tau - (k-1)*Tper.
%   5. Mark flyback samples NaN:  idx(ph >= T) = NaN  when Tfly > 0.
%
% INPUTS
%   t_abs  (double array) — absolute time samples [s]
%   t0     (double)       — CPI / snapshot start time [s]
%   T      (double)       — active chirp (ramp) duration [s]
%   Tfly   (double)       — flyback duration [s]; 0 for no flyback
%
% OUTPUT
%   idx    (double array) — chirp index array, same size as t_abs;
%                           NaN where t_abs falls in the flyback interval

    Tper = T + Tfly;
    tau  = t_abs - t0;

    k   = floor(tau ./ Tper);     % zero-based chirp count
    ph  = tau - k .* Tper;        % phase within the current PRI [s]
    idx = k + 1;                  % convert to 1-based chirp index

    % Mark flyback samples as NaN to exclude them from chirp masks.
    if Tfly > 0
        idx(ph >= T) = NaN;
    end
end