function ctx = ARIFOMA_logCpiInterferenceMetrics( ...
    ctx, idxTime, idxCpi, t0_s, ...
    intfDurAll_s, numInterferedPri, ratioPriInterfered, ...
    intfDurLos_s, intfDurRefl_s, intfRangeLos_m, intfRangeRefl_m, ...
    numInterferersPerPri, numLosInterferersPerPri, numReflInterferersPerPri, ...
    intfCenterFreq_Hz)
%LOGCPIINTERFERENCEMETRICS  Store interference metrics for one CPI in the context log.
%
% This function is the single point of truth for writing per-CPI metrics into ctx.
% It is intentionally strict: required fields/slots must already exist (preallocated
% during initialization). If the log schema is broken, this should error immediately.
%
% Inputs (per CPI):
%   idxTime                    Time-step index in the scenario/simulation.
%   idxCpi                     CPI index within the current time-step (micro-step).
%   t0_s                       Start time of the current CPI [s].
%
%   intfDurAll_s               Vector of interference event durations within the CPI [s].
%   numInterferedPri           Number of PRIs affected by interference within the CPI [-].
%   ratioPriInterfered         Fraction of PRIs interfered within the CPI (0..1) [-].
%
%   intfDurLos_s               Vector of durations for LOS-caused interference events [s].
%   intfDurRefl_s              Vector of durations for reflected-path interference events [s].
%
%   intfRangeLos_m             Vector of LOS path lengths associated with events [m].
%   intfRangeRefl_m            Vector of reflected path lengths associated with events [m].
%
%   numInterferersPerPri       Vector: number of interferers present per PRI [-].
%   numLosInterferersPerPri    Vector: LOS-only interferers per PRI [-].
%   numReflInterferersPerPri   Vector: reflected-only interferers per PRI [-].
%
%   intfCenterFreq_Hz          Optional scalar descriptor: center frequency of detected
%                              interference (e.g., "center-of-interference") [Hz].
%
% Logged schema (ctx.dataLog(idxTime).frame(idxCpi)):
%   t0_s
%   numInterferedPri
%   ratioPriInterfered
%   intfDurAll_s
%   intfDurLos_s
%   intfDurRefl_s
%   intfRangeLos_m
%   intfRangeRefl_m
%   numInterferersPerPri
%   numLosInterferersPerPri
%   numReflInterferersPerPri
%   intfCenterFreq_Hz
%   numLosEvents
%   numReflEvents

    % Store CPI reference time (used for alignment across time-steps and plotting).
    ctx.t0_s(1, idxTime) = t0_s;

    % Create a fresh per-CPI record (template guarantees consistent field presence).
    F = cpiMetricsTemplate();

    % PRI-level interference summary (dimensionless KPI for the CPI).
    F.numInterferedPri   = numInterferedPri;
    F.ratioPriInterfered = ratioPriInterfered;

    % Event-level interference durations (full set and split by propagation mechanism).
    F.intfDurAll_s  = intfDurAll_s;
    F.intfDurLos_s  = intfDurLos_s;
    F.intfDurRefl_s = intfDurRefl_s;

    % Corresponding path lengths (geometry-derived) for LOS and reflected components.
    F.intfRangeLos_m  = intfRangeLos_m;
    F.intfRangeRefl_m = intfRangeRefl_m;

    % Instantaneous interferer counts (per PRI), including LOS/reflection partitions.
    F.numInterferersPerPri     = numInterferersPerPri;
    F.numLosInterferersPerPri  = numLosInterferersPerPri;
    F.numReflInterferersPerPri = numReflInterferersPerPri;

    % Optional scalar descriptor of where interference energy is concentrated in frequency.
    F.intfCenterFreq_Hz = intfCenterFreq_Hz;

    % Compact scalar summaries (useful for plots/tables in the paper).
    F.numLosEvents  = numel(intfDurLos_s);
    F.numReflEvents = numel(intfDurRefl_s);

    % Commit to the preallocated log slot.
    ctx.dataLog(idxTime).frame(idxCpi) = F;
end


function F = cpiMetricsTemplate()
%CPI metrics template for strict logging (fixed schema).
%
% Notes:
% - Scalars default to NaN to make missing values obvious in analysis.
% - Vectors default to [].

F = struct( ...
    ... % PRI-level KPIs
    'numInterferedPri',           NaN, ...
    'ratioPriInterfered',         NaN, ...
    ... % Event-level durations (s)
    'intfDurAll_s',               [], ...
    'intfDurLos_s',               [], ...
    'intfDurRefl_s',              [], ...
    ... % Event-level path lengths (m)
    'intfRangeLos_m',             [], ...
    'intfRangeRefl_m',            [], ...
    ... % Interferer counts per PRI (-)
    'numInterferersPerPri',       [], ...
    'numLosInterferersPerPri',    [], ...
    'numReflInterferersPerPri',   [], ...
    ... % Frequency descriptor (Hz)
    'intfCenterFreq_Hz',          NaN, ...
    ... % Scalar event counts (-)
    'numLosEvents',               0, ...
    'numReflEvents',              0 );
end