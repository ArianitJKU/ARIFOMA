function Param = ARIFOMA_build_parameter_table(vehicleIds, Scenario, cfg)
% ARIFOMA_BUILD_PARAMETER_TABLE
%   Create the per-sensor radar parameter table for all vehicles in the scene.
%
% PURPOSE
% -------
%   This function is the single entry point for parameter generation in
%   ARIFOMA.  It performs two responsibilities:
%
%     1. Enriches cfg with scenario-derived fields that depend on the
%        loaded traffic data (e.g. per-vehicle corner inversion based on
%        lane position).
%
%     2. Dispatches to the appropriate parameter generator based on the
%        mode selected in cfg.paramMode.
%
%   By centralising the dispatch here, the rest of the simulation pipeline
%   (main_sim.m, ARIFOMA_prepare_scene_views, etc.) remains agnostic to
%   which generator is active.
%
% INPUTS
%   vehicleIds  (double array) — vector of unique vehicle IDs present in
%                                the scenario.
%   Scenario    (struct)       — loaded scenario struct; must contain at
%                                least Scenario.vLoc_table with columns
%                                Vid and V_y.
%   cfg         (struct)       — global configuration struct; key fields:
%
%       cfg.paramMode    (string) — parameter-generation mode; one of:
%                                    "EXTENDED"            per-corner random
%                                    "COMPASS"             direction-aware
%                                    "COMPASS_PORTIONS"    compass + partial band
%                                    "FIXED PARAMETER POOL" deterministic pool
%       cfg.lane_edges   (double array) — lane boundary positions [m],
%                                         used to assign vehicles to lanes
%                                         and determine corner inversion.
%
%   All other cfg fields required by the dispatched generator are passed
%   through unchanged.
%
% OUTPUT
%   Param  (table) — parameter table with one row per (vehicle, corner)
%                    sensor instance.  Schema matches the output of the
%                    dispatched generator (see individual generator
%                    functions for full column documentation).
%
% PARAMETER MODES
% ---------------
%   "EXTENDED"
%       Calls ARIFOMA_generate_corner_parameters.
%       Each corner draws waveform parameters independently and randomly
%       from the configured RF band and timing distributions.
%
%   "COMPASS"
%       Calls ARIFOMA_generate_corner_parameters_compass.
%       Carrier frequencies are assigned by the compass-facing direction
%       of each sensor corner, partitioning the 77-81 GHz band into four
%       1 GHz sub-bands (NW / NE / SE / SW).
%
%   "COMPASS_PORTIONS"
%       Calls ARIFOMA_generate_corner_parameters_compass.
%       As COMPASS but with configurable fractional band allocation per
%       sector (cfg.compassShare).
%
%   "FIXED PARAMETER POOL"
%       Calls ARIFOMA_generate_fixed_parameter_pool.
%       All vehicles draw from a pre-generated pool of exactly
%       cfg.numFixedWaveforms distinct configurations.  Enabled
%       automatically when cfg.UseFixedSeed = true.
%
% NOTES
% -----
%   - cfg is modified locally (cfg.CornerInvertMap is added) but the
%     caller's copy is not affected — MATLAB passes structs by value.
%   - Lane inversion: vehicles in lanes 4-6 have their left/right corner
%     labels swapped to reflect the opposing traffic direction.
%   - Any cfg enrichment that depends on Scenario should be added in this
%     function, not inside the individual generators.
%
% SEE ALSO
%   ARIFOMA_generate_corner_parameters,
%   ARIFOMA_generate_corner_parameters_compass,
%   ARIFOMA_generate_fixed_parameter_pool,
%   main_sim

% =========================================================================

%% ── Step 1: Scenario-Derived Configuration Enrichment ────────────────────
%   Compute per-vehicle corner inversion flags from lane position.
%   This must happen here because it requires Scenario.vLoc_table, which
%   is not available inside the individual generator functions.
%
%   cfg.CornerInvertMap is a table with columns {VehID, Invert} that the
%   compass generator uses to swap FL↔FR and RL↔RR for oncoming vehicles.

cfg.CornerInvertMap = make_corner_invert_map( ...
    Scenario.vLoc_table, cfg.lane_edges);

%% ── Step 2: Parameter Generator Dispatch ────────────────────────────────
%   Route to the appropriate generator based on cfg.paramMode.
%   mode is upper-cased to make matching case-insensitive.

mode = upper(string(cfg.paramMode));

switch mode

    case "EXTENDED"
        % Per-corner independent random parameter generation.
        % No compass or directional band constraints applied.
        Param = ARIFOMA_generate_corner_parameters(vehicleIds, cfg);

    case {"COMPASS", "COMPASS_PORTIONS"}
        % Direction-aware carrier frequency assignment.
        % COMPASS_PORTIONS additionally uses cfg.compassShare for partial
        % band allocation per sector.
        Param = ARIFOMA_generate_corner_parameters_compass(vehicleIds, cfg);

    case "FIXED PARAMETER POOL"
        % Constrained deployment: all vehicles draw from a pre-generated
        % pool of exactly cfg.numFixedWaveforms distinct configurations.
        Param = ARIFOMA_generate_fixed_parameter_pool(vehicleIds, cfg);

    otherwise
        error('ARIFOMA:buildParameterTable:unknownMode', ...
              ['Unknown cfg.paramMode = "%s". ' ...
               'Valid options: "EXTENDED", "COMPASS", ' ...
               '"COMPASS_PORTIONS", "FIXED PARAMETER POOL".'], ...
              cfg.paramMode);
end

end   % ARIFOMA_build_parameter_table

% =========================================================================
%% ── Local Helper Functions ───────────────────────────────────────────────
% =========================================================================

% -------------------------------------------------------------------------
function CornerInvertMap = make_corner_invert_map(vLoc_table, lane_edges)
% MAKE_CORNER_INVERT_MAP
%   Determine which vehicles require corner-label inversion based on their
%   lateral (Y) position in the traffic scene.
%
% BACKGROUND
% ----------
%   On a bidirectional highway, vehicles in opposing lanes travel in the
%   opposite direction.  From the perspective of a victim radar, the
%   "front-left" corner of an oncoming vehicle corresponds to the
%   "front-right" corner in the global coordinate frame.  Lanes 4-6
%   (the upper half of the lane_edges array) are treated as the opposing
%   direction and receive an inversion flag.
%
% ALGORITHM
% ---------
%   1. Compute the median lateral position of each vehicle across all
%      recorded time steps (robust to momentary lane changes).
%   2. Assign each vehicle to its nearest lane by minimising the distance
%      to the entries in lane_edges.
%   3. Set Invert = true for vehicles assigned to lanes 4, 5, or 6.
%
% INPUTS
%   vLoc_table  (table)        — vehicle location table with columns:
%                                  Vid  (vehicle ID)
%                                  V_y  (lateral position [m])
%   lane_edges  (double array) — lateral positions of lane centres [m],
%                                e.g. [4 8 12 16 20 24] for a 6-lane road.
%
% OUTPUT
%   CornerInvertMap  (table) — two-column table with:
%                                VehID  (double)  — vehicle ID
%                                Invert (logical) — true if corners are inverted

    vehIds = unique(vLoc_table.Vid(:));
    nV     = numel(vehIds);

    laneIdx = nan(nV, 1);
    y_med   = nan(nV, 1);

    % ── Assign each vehicle to its nearest lane ───────────────────────────
    for i = 1:nV
        rows      = (vLoc_table.Vid == vehIds(i));
        y_med(i)  = median(vLoc_table.V_y(rows), 'omitnan');

        % Nearest-neighbour assignment to lane_edges.
        [~, laneIdx(i)] = min(abs(lane_edges - y_med(i)));
    end

    % ── Flag opposing-direction lanes (4, 5, 6) for corner inversion ──────
    %   Lanes 1-3 are treated as the ego direction (no inversion).
    %   Lanes 4-6 are treated as the opposing direction (invert corners).
    Invert = ismember(laneIdx, [4 5 6]);

    CornerInvertMap = table(vehIds, Invert, ...
                            'VariableNames', {'VehID', 'Invert'});

end   % make_corner_invert_map