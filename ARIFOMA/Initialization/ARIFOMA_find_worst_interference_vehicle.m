function worst = ARIFOMA_find_worst_interference_vehicle( ...
        vLoc, intf, Param, cfg, opt)
% ARIFOMA_FIND_WORST_INTERFERENCE_VEHICLE
%   Identify the sensor-equipped vehicle that experiences the highest
%   interference load across the simulation time window.
%
% PURPOSE
% -------
%   Scans the pre-computed interference path table (intf) over the
%   acquisition time grid defined by cfg.simTime and ranks all
%   sensor-equipped victim vehicles by their interference load.  The
%   vehicle with the highest load is selected as the worst-case victim and
%   a full scene frame is extracted at its peak-interference time.
%
%   The result is used by main_sim.m when cfg.worstVehicle = true to
%   automatically select the most informative victim for analysis, rather
%   than relying on a manually specified or randomly drawn vehicle ID.
%
% RANKING METRICS
% ---------------
%   Two metrics are supported, selected via cfg.worstMetric:
%
%     "peak"  (default) — vehicle with the highest number of simultaneous
%                         interference paths in any single acquisition frame.
%                         Reflects the worst-case instantaneous interference
%                         density experienced by that vehicle.
%
%     "total" — vehicle with the highest cumulative interference path count
%               across the entire time window.  Reflects sustained exposure
%               rather than peak load.
%
% ALGORITHM
% ---------
%   1. Validate required columns in intf and Param.
%   2. Build the acquisition time grid from cfg.simTime.
%   3. Restrict analysis to vehicles with HasSensor = true in Param.
%   4. Filter intf rows to the time window and to eligible victim vehicles.
%   5. For each eligible victim vehicle, compute:
%        TotalPaths — total interference path count across all frames
%        PeakPaths  — maximum path count in any single frame
%        tPeak      — acquisition time at which PeakPaths occurred
%   6. Select the worst vehicle by the chosen metric.
%   7. Extract a scene frame at tWorst using ARIFOMA_frame_extract_mc.
%   8. Pack results into the worst output struct.
%
% INPUTS
%   vLoc  (table)  — vehicle location table (Scenario.vLoc_table)
%   intf  (table)  — interference path table (Scenario.intf_table);
%                    required columns: time, interferedId, interferingId
%   Param (table)  — ARIFOMA parameter table; required columns:
%                    VehID, HasSensor
%   cfg   (struct) — global configuration struct; required fields:
%                      cfg.simTime      (double array) — acquisition time
%                                        grid [s]; used as the ranking
%                                        time window
%                    optional fields:
%                      cfg.worstMetric  (string) — "peak" (default) or
%                                        "total"
%   opt   (struct) — options struct passed through to
%                    ARIFOMA_frame_extract_mc for frame extraction
%
% OUTPUT
%   worst  (struct) — worst-case vehicle summary with fields:
%     .vehIdWorst       (double) — vehicle ID of the worst-case victim
%     .tWorst           (double) — acquisition time of peak interference [s]
%     .metric           (char)   — metric used: 'peak' or 'total'
%     .worstMetricValue (double) — metric value for the worst vehicle
%     .tablePerVeh      (table)  — per-vehicle ranking table, sorted
%                                  descending by the chosen metric;
%                                  columns: VehId, TotalPaths, PeakPaths,
%                                  tPeak
%     .frame            (struct) — scene frame struct extracted at tWorst
%                                  for vehIdWorst (output of
%                                  ARIFOMA_frame_extract_mc)
%
% NOTES
% -----
%   - Only vehicles with HasSensor = true in Param are eligible as victims.
%     Vehicles without sensors cannot experience radar interference.
%   - Interference rows in intf are binned to the nearest acquisition time
%     index using round((t - tStart) / dt).  This assumes intf.time aligns
%     approximately with cfg.simTime; small misalignments are absorbed by
%     the rounding.
%   - The median step size of cfg.simTime is used as dt.  If the grid has
%     fewer than 2 points or produces a non-positive dt, dt defaults to 0.1 s.
%
% SEE ALSO
%   ARIFOMA_frame_extract_mc,
%   ARIFOMA_build_parameter_table,
%   main_sim
 
% =========================================================================
 
%% ── Step 1: Input Validation ─────────────────────────────────────────────
%   Fail early with informative errors if required table columns are absent.
 
mustHaveCols(intf, {'time','interferedId','interferingId'});
 
if ~ismember('VehID',     Param.Properties.VariableNames) || ...
   ~ismember('HasSensor', Param.Properties.VariableNames)
    error('ARIFOMA:findWorstVehicle:missingColumn', ...
          'Param must contain columns VehID and HasSensor.');
end
 
%% ── Step 2: Build Acquisition Time Grid ──────────────────────────────────
%   Use cfg.simTime as the ranking time window.  This matches the main
%   simulation grid so results are directly comparable to a full run.
 
if ~isfield(cfg, 'simTime') || isempty(cfg.simTime)
    error('ARIFOMA:findWorstVehicle:noTimeGrid', ...
          'cfg.simTime is required and must be non-empty.');
end
 
tGrid  = double(cfg.simTime(:));
tStart = tGrid(1);
tEnd   = tGrid(end);
nT     = numel(tGrid);
 
% Derive step size robustly; default to 0.1 s if grid has fewer than 2 points.
dt = median(diff(tGrid));
if ~isfinite(dt) || dt <= 0
    dt = 0.1;
end
 
%% ── Step 3: Identify Sensor-Equipped Vehicles ────────────────────────────
%   Only vehicles with at least one active sensor corner are eligible as
%   victims.  Vehicles without sensors cannot receive interference.
 
vehWithSensors = unique(double(Param.VehID(Param.HasSensor)));
vehWithSensors = vehWithSensors(isfinite(vehWithSensors));
 
if isempty(vehWithSensors)
    error('ARIFOMA:findWorstVehicle:noSensorVehicles', ...
          'No sensor-equipped vehicles found in Param.');
end
 
%% ── Step 4: Filter Interference Table ───────────────────────────────────
%   Retain only rows within the time window and where the victim vehicle
%   carries at least one sensor.
 
ti    = double(intf.time);
inWin = (ti >= tStart) & (ti <= tEnd);
Rwin  = intf(inWin, :);
tiw   = double(Rwin.time);
 
if isempty(Rwin)
    error('ARIFOMA:findWorstVehicle:noRowsInWindow', ...
          'No interference rows in time window [%.3f, %.3f] s.', tStart, tEnd);
end
 
% Map each row's timestamp to the nearest time-grid index (1-based).
k = round((tiw - tStart) / dt) + 1;
k = max(1, min(nT, k));
 
% Keep only rows where the victim is sensor-equipped.
victimIds        = double(Rwin.interferedId);
isVictimEligible = ismember(victimIds, vehWithSensors);
Rwin             = Rwin(isVictimEligible, :);
k                = k(isVictimEligible);
 
if isempty(Rwin)
    error('ARIFOMA:findWorstVehicle:noEligibleRows', ...
          'No interference rows found for sensor-equipped vehicles in the selected window.');
end
 
%% ── Step 5: Per-Vehicle Interference Metrics ─────────────────────────────
%   For each eligible victim vehicle, compute TotalPaths, PeakPaths, and
%   the acquisition time at which the peak occurred.
 
vehUnique = unique(double(Rwin.interferedId));
vehUnique = vehUnique(isfinite(vehUnique));
nVeh      = numel(vehUnique);
 
totalCount = zeros(nVeh, 1);
peakCount  = zeros(nVeh, 1);
tAtPeak    = nan(nVeh, 1);
 
intfVict = double(Rwin.interferedId);
 
for ii = 1:nVeh
    vid   = vehUnique(ii);
    rowsV = (intfVict == vid);
 
    % Total interference path count across the full time window.
    totalCount(ii) = sum(rowsV);
 
    % Per-frame path count via accumarray; peak = worst single frame.
    kk        = k(rowsV);
    cPerFrame = accumarray(kk, 1, [nT 1], @sum, 0);
 
    [pk, idxPk]  = max(cPerFrame);
    peakCount(ii) = pk;
    tAtPeak(ii)   = tGrid(idxPk);
end
 
Tveh = table(vehUnique, totalCount, peakCount, tAtPeak, ...
    'VariableNames', {'VehId','TotalPaths','PeakPaths','tPeak'});
 
%% ── Step 6: Select Worst Vehicle by Chosen Metric ────────────────────────
%   Default metric is "peak".  cfg.worstMetric overrides if present.
 
metric = "peak";
if isfield(cfg, 'worstMetric') && ~isempty(cfg.worstMetric)
    metric = lower(string(cfg.worstMetric));
end
 
switch metric
    case "total"
        [~, ixW]   = max(Tveh.TotalPaths);
        worstValue = Tveh.TotalPaths(ixW);
    otherwise
        % Default: peak simultaneous path count.
        [~, ixW]   = max(Tveh.PeakPaths);
        worstValue = Tveh.PeakPaths(ixW);
        metric     = "peak";
end
 
vehIdWorst = Tveh.VehId(ixW);
tWorst     = Tveh.tPeak(ixW);
 
fprintf('[ARIFOMA] Worst vehicle: ID=%g | metric=%s | value=%g | tWorst=%.3f s\n', ...
    vehIdWorst, metric, worstValue, tWorst);
 
%% ── Step 7: Extract Scene Frame at Worst Time ────────────────────────────
%   Extract the full scene frame at tWorst for the worst-case vehicle so
%   that main_sim.m can initialise its simulation from this starting point.
 
opt2           = opt;
opt2.VehicleId = vehIdWorst;   % focus frame extraction on the worst victim
frame          = ARIFOMA_frame_extract_mc(vLoc, intf, Param, tWorst, opt2);
 
%% ── Step 8: Pack Output Struct ───────────────────────────────────────────
worst                  = struct();
worst.vehIdWorst       = vehIdWorst;
worst.tWorst           = tWorst;
worst.metric           = char(metric);
worst.worstMetricValue = worstValue;
worst.tablePerVeh      = sortrows(Tveh, metricSortVar(metric), 'descend');
worst.frame            = frame;
 
end   % ARIFOMA_find_worst_interference_vehicle
 
% =========================================================================
%% ── Local Helper Functions ───────────────────────────────────────────────
% =========================================================================
 
% -------------------------------------------------------------------------
function vname = metricSortVar(metric)
% METRICSORTVAR  Return the Tveh column name corresponding to a metric string.
%
%   Used to select the sort variable when producing worst.tablePerVeh.
%
%   INPUT
%     metric  (string) — "total" or "peak"
%
%   OUTPUT
%     vname   (char)   — column name: 'TotalPaths' or 'PeakPaths'
 
    if metric == "total"
        vname = 'TotalPaths';
    else
        vname = 'PeakPaths';
    end
end
 
% -------------------------------------------------------------------------
function mustHaveCols(T, cols)
% MUSTHAVECOLS  Assert that a table contains all required column names.
%
%   Throws an informative error immediately if any column is absent,
%   preventing downstream indexing errors that are harder to diagnose.
%
%   INPUTS
%     T     (table)     — table to validate
%     cols  (cell array of char) — required column names
 
    for i = 1:numel(cols)
        if ~ismember(cols{i}, T.Properties.VariableNames)
            error('ARIFOMA:mustHaveCols:missingColumn', ...
                  'Required column "%s" is absent from the interference table.', ...
                  cols{i});
        end
    end
end
 