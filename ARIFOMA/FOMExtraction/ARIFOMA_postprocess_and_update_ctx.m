function ctx = ARIFOMA_postprocess_and_update_ctx( ...
        ctx, out, corner, cfg, kCPI, nCPI, kT, t0_s)
% ARIFOMA_POSTPROCESS_AND_UPDATE_CTX
%   Finalise one CPI micro-step and log all interference figures of merit.
%
% PURPOSE
% -------
%   This function is the final stage of the per-corner CPI processing
%   pipeline.  It transfers the updated context from the strategy output,
%   maps the raw metric arrays returned by the strategy into the labelled
%   FOM struct defined in ARIFOMA, commits the FOM struct to the persistent
%   data log, and triggers end-of-CPI updates for any adaptive scheme logic
%   (BATS victim and interferer statistics).
%
% PIPELINE POSITION
% -----------------
%   Called from runCPI_perCorner (main_sim.m) as Stage 6 (final stage):
%
%     Stage 5                              Stage 6
%   ARIFOMA_run_strategy  →  ARIFOMA_postprocess_and_update_ctx
%
% FIGURES OF MERIT (FOMs)
% -----------------------
%   The FOM struct F logged to ctx.dataLog at each (kT, kCPI) pair maps
%   directly to the unified interference-related FOMs defined in Table 1
%   of the ARIFOMA paper:
%
%   FOM1 — Per-CPI interference ratio
%     F.FOM1_numInterferedPri     number of victim PRIs with any overlap
%     F.FOM1_ratioInterferedPri   fraction of victim PRIs interfered [0,1]
%
%   FOM2 — Normalised interference duration  (T_I / T_chirp per occurrence)
%     F.FOM2_intfDurAll_norm      combined LOS + reflected durations
%     F.FOM2_intfDurLos_norm      LOS durations only
%     F.FOM2_intfDurRefl_norm     reflected durations only
%
%   FOM3 — Interference path origin and distance
%     F.FOM3_intfRangeLos_m       LOS path distances [m]
%     F.FOM3_intfRangeRefl_m      total reflected path distances [m]
%     F.FOM3_intfRangeRefl_d1_m   first reflected leg distances [m]
%     F.FOM3_intfRangeRefl_d2_m   second reflected leg distances [m]
%     F.FOM3_numLosEvents         number of LOS interference occurrences
%     F.FOM3_numReflEvents        number of reflected interference occurrences
%
%   FOM4 — Number of interferers per PRI
%     F.FOM4_numIntfPerPri        total interferers per victim PRI
%     F.FOM4_numLosIntfPerPri     LOS interferers per victim PRI
%     F.FOM4_numReflIntfPerPri    reflected interferers per victim PRI
%
%   BATS observable
%     F.centerofIntf              duration-weighted interference centroid
%                                 estimate r_hat ∈ [0,1]; NaN if no
%                                 interference or BATS not active
%
% INPUTS
%   ctx     (struct)  — current corner context; modified in-place
%   out     (struct)  — strategy output struct from ARIFOMA_run_strategy;
%                       required fields:
%                         out.ctx     — updated context from strategy
%                         out.metrics — per-CPI metric arrays (M)
%   corner  (struct)  — corner bundle from ARIFOMA_apply_all_schemes;
%                       required fields:
%                         corner.V.Scheme         — victim scheme string
%                         corner.schemeInterferer  — interferer scheme string
%                         corner.I                 — interferer struct array
%   cfg     (struct)  — global configuration struct; retained for interface
%                       compatibility but not used directly in this function
%   kCPI    (int)     — current CPI micro-step index (1-based)
%   nCPI    (int)     — total CPI micro-steps per acquisition snapshot
%   kT      (int)     — current acquisition time index (1-based)
%   t0_s    (double)  — current acquisition timestamp [s]
%
% OUTPUT
%   ctx  (struct) — updated corner context with:
%                     ctx.t0(1,kT)             — acquisition timestamp [s]
%                     ctx.frameHadInterfV       — accumulated interference
%                                                flag for current snapshot
%                     ctx.wasInterferedPrevFrameV — interference flag from
%                                                   completed snapshot
%                     ctx.dataLog(kT).CPI(kCPI) — logged FOM struct F
%
% NOTES
% -----
%   - cfg is accepted as an argument for interface compatibility with the
%     runCPI_perCorner call signature but is not used within this function.
%   - ctx.frameHadInterfV accumulates across CPI micro-steps within a
%     single acquisition snapshot.  It is committed to
%     ctx.wasInterferedPrevFrameV only at the last CPI (kCPI == nCPI) so
%     that adaptive schemes can query whether the previous complete snapshot
%     contained any interference.
%   - The BATS end-of-CPI updates are triggered unconditionally; the
%     called functions internally check whether BATS is active for the
%     victim and interferer respectively.
%
% SEE ALSO
%   ARIFOMA_run_strategy,
%   ARIFOMA_update_victim_bats_stats_end_cpi,
%   ARIFOMA_update_interferer_bats_stats_end_cpi,
%   runCPI_perCorner (main_sim.m)

% =========================================================================

%% ── Transfer Updated Context from Strategy Output ────────────────────────
%   The strategy (ARIFOMA_run_strategy) may modify ctx internally.
%   The updated copy is returned in out.ctx and must be adopted here before
%   any further context writes.

ctx = out.ctx;

%% ── Extract Per-CPI Metric Arrays ───────────────────────────────────────
%   M is the raw metric struct returned by the active strategy.
%   All FOM fields are derived from M in the logging section below.

M = out.metrics;

%% ── Interference State Accumulation ─────────────────────────────────────
%   hadIntf flags whether the current CPI micro-step contains any
%   interference (at least one non-zero normalised duration).
%
%   ctx.frameHadInterfV accumulates this flag across all CPI micro-steps
%   within the current acquisition snapshot using logical OR.  It is
%   committed to ctx.wasInterferedPrevFrameV at the final micro-step
%   (kCPI == nCPI) so that adaptive schemes can query the previous
%   snapshot's interference status at the start of the next snapshot.

hadIntf = any(M.all_interference_durations > 0);
ctx.frameHadInterfV = ctx.frameHadInterfV || hadIntf;

if kCPI == nCPI
    % Final CPI micro-step: commit the snapshot-level interference flag.
    ctx.wasInterferedPrevFrameV = ctx.frameHadInterfV;
end

%% ── Store Acquisition Timestamp ──────────────────────────────────────────
ctx.t0(1, kT) = t0_s;

%% ── Build FOM Log Entry ──────────────────────────────────────────────────
%   Initialise the FOM struct with NaN / empty sentinels so that partially
%   populated entries (e.g. when no interference is detected) still have
%   a consistent schema in the data log.

F = struct( ...
    'FOM1_numInterferedPri',      NaN, ...
    'FOM1_ratioInterferedPri',    NaN, ...
    'FOM2_intfDurAll_norm',       [], ...
    'FOM2_intfDurLos_norm',       [], ...
    'FOM2_intfDurRefl_norm',      [], ...
    'FOM3_intfRangeLos_m',        [], ...
    'FOM3_intfRangeRefl_m',       [], ...
    'FOM3_intfRangeRefl_d1_m',    [], ...
    'FOM3_intfRangeRefl_d2_m',    [], ...
    'FOM3_numLosEvents',          0, ...
    'FOM3_numReflEvents',         0, ...
    'FOM4_numIntfPerPri',         [], ...
    'FOM4_numLosIntfPerPri',      [], ...
    'FOM4_numReflIntfPerPri',     [], ...
    'centerofIntf',               NaN);

%% ── Assign FOM Values from Strategy Metrics ──────────────────────────────
%   Map the raw metric arrays in M to the labelled FOM fields.
%   FOM numbering follows Table 1 of the ARIFOMA paper.

% FOM1 — Per-CPI interference ratio.
F.FOM1_numInterferedPri   = M.total_interfered_chirps;
F.FOM1_ratioInterferedPri = M.P_coll_per_n;

% FOM2 — Normalised interference duration (T_I / T_chirp per occurrence).
F.FOM2_intfDurAll_norm   = M.all_interference_durations;
F.FOM2_intfDurLos_norm   = M.all_interference_durations_LOS;
F.FOM2_intfDurRefl_norm  = M.all_interference_durations_REFL;

% FOM3 — Interference path origin and distance.
F.FOM3_intfRangeLos_m     = M.all_interference_dist_LOS;
F.FOM3_intfRangeRefl_m    = M.all_interference_dist_REFL;
F.FOM3_intfRangeRefl_d1_m = M.all_interference_dist_REFL_d1;
F.FOM3_intfRangeRefl_d2_m = M.all_interference_dist_REFL_d2;
F.FOM3_numLosEvents        = numel(M.all_interference_durations_LOS);
F.FOM3_numReflEvents       = numel(M.all_interference_durations_REFL);

% FOM4 — Number of distinct interferers per victim PRI.
F.FOM4_numIntfPerPri      = M.all_nIntf_per_chirp_total;
F.FOM4_numLosIntfPerPri   = M.all_nIntf_per_chirp_LOS;
F.FOM4_numReflIntfPerPri  = M.all_nIntf_per_chirp_REFL;

% BATS observable — interference centroid estimate r_hat.
% NaN when no interference is detected or BATS scheme is not active.
F.centerofIntf = M.all_r_hat_sub;

%% ── Commit FOM Entry to Data Log ─────────────────────────────────────────
%   The log is indexed by acquisition time step (kT) and CPI micro-step
%   (kCPI), forming a 2D grid of FOM entries across the full simulation.

ctx.dataLog(kT).CPI(kCPI) = F;

%% ── End-of-CPI Adaptive Scheme Updates ───────────────────────────────────
%   Trigger BATS statistic updates for the victim and all interferers.
%   These functions internally check whether BATS is active for each sensor
%   before updating, so they are safe to call unconditionally.
%
%   Updates use the completed CPI's FOM log (now committed above) to
%   compute the updated interference centroid estimate for the next CPI.

ctx = ARIFOMA_update_victim_bats_stats_end_cpi( ...
    ctx, string(corner.V.Scheme), kT, kCPI, nCPI);

ctx = ARIFOMA_update_interferer_bats_stats_end_cpi( ...
    ctx, corner.schemeInterferer, corner.I, kCPI, nCPI);

end   % ARIFOMA_postprocess_and_update_ctx