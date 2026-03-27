%% main_sim.m — ARIFOMA
% =========================================================================
%
% Automotive Radar Interference Figures of Merit Analysis Framework
%
% PURPOSE
% -------
%   Simulates FMCW radar-level interference on a multi-lane highway
%   scene and evaluates interference mitigation strategies for
%   corner-mounted automotive radar sensors.
%
% SIMULATION PIPELINE
% -------------------
%   1. Load a pre-computed traffic / propagation scenario from disk.
%   2. Generate radar parameters for all sensor-equipped vehicles.
%   3. Select a victim vehicle (worst-case, specific ID, or random).
%   4. Iterate over acquisition time snapshots (cfg.simTime).
%      a. Extract and draw the per-snapshot scene geometry.
%      b. Build per-corner victim / interferer link tables.
%      c. For each CPI micro-step, run the selected waveform / mitigation
%         strategy and accumulate the resulting radar context.
%   5. Save the final `radarCtx` struct for post-processing and analysis.
%
% OUTPUT
% ------
%   radarCtx<cfg.sceneTag>.mat — per-corner radar context struct containing
%   all CPI logs, metrics, and the victim parameter snapshot.
%
% DEPENDENCIES
% ------------
%   loadScenario                    (local helper, this file)
%   ARIFOMA_build_parameter_table
%   ARIFOMA_find_worst_interference_vehicle
%   ARIFOMA_apply_victim_overrides
%   ARIFOMA_init_radarcorner_contexts
%   prepopulate_lastF0BySensor
%   ARIFOMA_frame_extract_mc
%   ARIFOMA_frame_draw_mc
%   ARIFOMA_prepare_scene_views
%   ARIFOMA_apply_all_schemes
%   ARIFOMA_build_engine_inputs
%   ARIFOMA_run_strategy
%   ARIFOMA_postprocess_and_update_ctx
%
% NOTES
% -----
%   • The simulation is built around per-corner victim / interferer views;
%     each vehicle corner (FL, FR, RL, RR) is treated independently.
%   • Parameter-generation mode is controlled by cfg.paramMode.
%     Supported modes: "COMPASS", "EXTENDED", "COMPASS_PORTIONS",
%     "FIXED PARAMETER POOL".
%   • All results accumulate inside `radarCtx`, which is saved at the end
%     and can be fed directly into ARIFOMA post-processing scripts.
%
% USAGE
%   Run as a MATLAB script.  Adjust cfg fields in the "User Configuration"
%   section before execution.
%
% AUTHORS
%   Arianit Preniqi, Oliver Lang, Stefan Schmalzl, Reinhard Feger
%
%   Johannes Kepler University Linz, Linz, Austria
% VERSION
%   v1.0  —  26.03.2026
%
% LICENSE
%   [License statement — fill in before publication]
%
% =========================================================================
 
%% ── Workspace Initialisation ─────────────────────────────────────────────
clear; close all; clc;
 
% ── Global figure defaults (LaTeX rendering, consistent font sizes) ───────
set(groot,'defaultTextInterpreter','latex', ...
          'defaultAxesTickLabelInterpreter','latex', ...
          'defaultLegendInterpreter','latex', ...
          'defaultColorbarTickLabelInterpreter','latex');
set(groot, ...
    'defaultAxesFontSize',14, ...
    'defaultTextFontSize',14, ...
    'defaultLegendFontSize',14, ...
    'defaultColorbarFontSize',14, ...
    'defaultAxesTitleFontSizeMultiplier',1, ...
    'defaultAxesLabelFontSizeMultiplier',1);
 
% =========================================================================
%% ── Section 1: User Configuration ───────────────────────────────────────
% =========================================================================
%   All experiment knobs are collected here so that no parameter is buried
%   inside the processing loop.  Sub-sections group related settings.
% =========================================================================
 
cfg = struct();
 
% ── 1.1  Simulation Time Axis ─────────────────────────────────────────────
%   Each entry in cfg.simTime defines one acquisition snapshot (seconds).
%   The victim scene is extracted and processed independently at each step.
cfg.TA=0.1;
cfg.TS=30;
cfg.simTime = 0:cfg.TA:cfg.TS;          % [s]  Acquisition times (0 – 30 s, 0.1 s steps)
 
% ── 1.2  Scenario / Output Naming ────────────────────────────────────────
%   cfg.sceneTag is appended to every output file produced by this run,
%   making results from different experiments easy to distinguish.
cfg.sceneTag    = '_SRR_25MHz_density85';
cfg.bats_eps     = 0.07;         %      BATS algorithm convergence threshold
cfg.outVideoFile = sprintf('scene%s.mp4', cfg.sceneTag);
 
% ── 1.3  Playback / Animation Settings ───────────────────────────────────
%   Used when video writing or frame hold is enabled.
cfg.fps            = 24;         % [fps] Output video frame rate
cfg.holdTopFrames  = 1;          % [frames] Hold count for top-panel frame
cfg.holdTileFrames = 2;          % [frames] Hold count for corner-tile frames
 
% ── 1.4  Chirp Flyback / Return-fly Timing Limits ────────────────────────
%   These bounds define legal flyback durations for waveform generation.
%   Stepped-flyback modes use tighter bounds than the default mode.
cfg.FlybackDefaultMin_s = 3e-6;  % [s]  Default flyback lower bound
cfg.FlybackDefaultMax_s = 10e-6; % [s]  Default flyback upper bound
cfg.FlybackSteppedMin_s = 2e-6;  % [s]  Stepped flyback lower bound
cfg.FlybackSteppedMax_s = 4e-6;  % [s]  Stepped flyback upper bound
cfg.ReturnFlyMin_s      = 5e-6;  % [s]  Return-fly lower bound
cfg.ReturnFlyMax_s      = 7e-6;  % [s]  Return-fly upper bound
 
% ── 1.5  Per-Class RF Band and Waveform Parameter Ranges ─────────────────
%   cfg.bandByScenario defines the operating frequency band, sweep
%   bandwidth range, and chirp period range for each radar class.
%   Classes: MRR (medium-range), SRR (short-range), LRR (long-range).
cfg.bandByScenario = struct();
 
cfg.bandByScenario.MRR = struct( ...
    'bandMin',           77e9, ...   % [Hz]  Lower RF bound
    'bandMax',           81e9, ...   % [Hz]  Upper RF bound
    'Brange',     [400e6 800e6], ... % [Hz]  Sweep bandwidth range
    'Trange',    [20e-6  40e-6], ... % [s]   Chirp period range
    'T_flyback_default', [3e-6 5e-6], ... % [s]  Default flyback range
    'T_flyback_stepped', [1e-6 3e-6], ... % [s]  Stepped flyback range
    'T_returnfly',       [4e-6 6e-6]);    % [s]  Return-fly range
 
cfg.bandByScenario.SRR = struct( ...
    'bandMin',            77e9, ...
    'bandMax',            81e9, ...
    'Brange',      [1e9    2e9], ...
    'Trange',    [10e-6  20e-6], ...
    'T_flyback_default', [3e-6 5e-6], ...
    'T_flyback_stepped', [1e-6 3e-6], ...
    'T_returnfly',       [4e-6 6e-6]);
 
cfg.bandByScenario.LRR = struct( ...
    'bandMin',            76e9, ...
    'bandMax',            77e9, ...
    'Brange',    [200e6 500e6], ...
    'Trange',    [40e-6 100e-6], ...
    'T_flyback_default', [3e-6 5e-6], ...
    'T_flyback_stepped', [1e-6 3e-6], ...
    'T_returnfly',       [4e-6 6e-6]);
 
% ── 1.6  Parameter-Generator Mode ────────────────────────────────────────
%   Controls how radar parameters are assigned to vehicles.
%   Recognised values:
%     "COMPASS"             — direction-aware frequency partitioning
%     "EXTENDED"            — extended random parameter space
%     "FIXED PARAMETER POOL"— deterministic pool (requires UseFixedSeed)
cfg.paramMode    = "EXTENDED";
cfg.UseFixedSeed = false;
cfg.numFixedWaveforms=10;
cfg.Seed        = 7;             %      RNG seed (active when UseFixedSeed = true)
cfg.pHasSensor  = 1;             %      Probability that a vehicle carries a sensor
 
% ── 1.7  Victim Corner Scenario Assignment ────────────────────────────────
%   Specifies the radar class assigned to each corner of the victim vehicle.
%   Used by the victim override pipeline.
cfg.victimscenarioPerCorner = struct( ...
    'FL', "MRR", ...
    'FR', "MRR", ...
    'RL', "MRR", ...
    'RR', "MRR");
 
% ── 1.8  Radar-Class Distribution for the General Traffic Population ──────
%   Per-corner probability weights for MRR / SRR / LRR assignment.
%   Weights do not need to sum to 1 (they are normalised internally).
cfg.scenarioDist.FL.list    = ["MRR","SRR","LRR"];
cfg.scenarioDist.FL.weights = [1 0 0];
 
cfg.scenarioDist.FR.list    = ["MRR","SRR","LRR"];
cfg.scenarioDist.FR.weights = [1 0 0];
 
cfg.scenarioDist.RL.list    = ["MRR","SRR","LRR"];
cfg.scenarioDist.RL.weights = [1 0 0];
 
cfg.scenarioDist.RR.list    = ["MRR","SRR","LRR"];
cfg.scenarioDist.RR.weights = [1 0 0];
 
% ── 1.9  Waveform and Scheme Pools for Traffic Population ─────────────────
%   Weights select which scheme / waveform / step parameters are drawn
%   during randomised population assignment.
cfg.SchemeList    = ["CONVENTIONAL","RANDOM_FH_PER_CPI","FAR", ...
                     "BLUEFMCW","BATSINSPIRED_FH_PER_CPI"];
cfg.SchemeWeights = [1 0 0 0 0];   % All traffic uses CONVENTIONAL
 
cfg.WaveformList    = ["POSITIVE-SLOPE","POSITIVE-SLOPE STEPPED", ...
                       "NEGATIVE-SLOPE","NEGATIVE-SLOPE STEPPED"];
cfg.WaveformWeights = [1 0 0 0];   % All traffic uses POSITIVE-SLOPE
 
cfg.StepCountMin = 64;             %      Minimum number of stepped-chirp sub-sweeps
cfg.StepCountMax = 128;            %      Maximum number of stepped-chirp sub-sweeps
cfg.StepHzMin    = 10e6;           % [Hz] Minimum frequency step size
cfg.StepHzMax    = 20e6;           % [Hz] Maximum frequency step size
 
% ── 1.10  ADC and Detector Settings ──────────────────────────────────────
cfg.BADC_Hz       = 25e6;          % [Hz] ADC sampling bandwidth
cfg.randstartTime = 40e-6;         % [s]  Maximum randomised chirp start offset
 
% ── 1.11  COMPASS Mode: Lane and Band Settings ────────────────────────────
%   COMPASS partitions the 77–81 GHz band across cardinal compass sectors
%   and assigns each sector to a lane group.
cfg.numberoflanes = 6;
cfg.lane_edges    = 4 + (0 : cfg.numberoflanes-1) * 4;   % [m] Lane boundary positions
cfg.compassShare  = 1;             %      Fraction of band used per compass sector
 
cfg.CPI_duration_s = 10e-3;       % [s]  Coherent processing interval (CPI) duration
 
% ── 1.12  Plotting and Debug Flags ───────────────────────────────────────
cfg.plotScenario                = false;  % Draw scenario at each snapshot
cfg.plotScenarioEnd             = true;   % Draw scenario at final snapshot only
cfg.boolIgnoreFlybackInterferer = true;   % Exclude flyback overlap from interference count
 
% ── 1.13  COMPASS Carrier-Frequency Band Boundaries ──────────────────────
%   Each cardinal direction maps to a 1 GHz sub-band within 77–81 GHz.
cfg.compassFcBands.NW = [77e9 78e9];   % [Hz]
cfg.compassFcBands.NE = [78e9 79e9];   % [Hz]
cfg.compassFcBands.SE = [79e9 80e9];   % [Hz]
cfg.compassFcBands.SW = [80e9 81e9];   % [Hz]
 
% ── 1.14  Victim Vehicle Selection ───────────────────────────────────────
%   Exactly one selection mode should be active at a time:
%     worstVehicle    — automatically find the highest-interference vehicle
%     specificVehicle — use the fixed ID in cfg.VehID
%     neither         — draw a random vehicle that carries at least one sensor
cfg.worstVehicle    = false;
cfg.specificVehicle = true;
cfg.VehID           = 243;         %      Victim vehicle ID (specificVehicle mode)
 
% ── 1.15  Victim Waveform / Scheme Override ───────────────────────────────
%   Applied by the victim override pipeline when the parameter mode permits.
cfg.adaptWaveformVictim = "POSITIVE-SLOPE";
cfg.adaptSchemeVictim   = "BATSINSPIRED_FH_PER_CPI";
 
% ── 1.16  CPI Micro-Steps per Acquisition Snapshot ───────────────────────
%   Each acquisition snapshot is split into this many CPI micro-steps,
%   allowing intra-snapshot adaptation to be captured.
cfg.numberofCPIsperAcquisition = 2;
 
% ── 1.17  Corner Ordering ─────────────────────────────────────────────────
%   Used consistently across the parameter table, link builder, and
%   context struct to index the four radar corners.
cfg.corners         = ["FL","RR","RL","FR"];
cfg.SensorIdMap     = table();     % Optional external sensor-ID mapping table
cfg.KeepIdForAbsent = true;        % Retain sensor IDs for corners with no active sensor
 
% ── 1.18  Input Scenario File ─────────────────────────────────────────────
%   Path to the pre-computed SUMO / WiLabV traffic and interference table.
cfg.inputMat = 'C:\Users\PC\Desktop\ARIFOMA\SUMO_WiLabVISIM_Data\Highway_Density85_3Lanes.mat';
 
% =========================================================================
%% ── Section 2: Figure and Tiled-Layout Initialisation ───────────────────
% =========================================================================
%   A 3×2 tiled layout is created once and reused across all time steps.
%   Top row (spanning both columns): global scene overview.
%   Lower 2×2 tiles (one per corner): per-corner victim / interferer view.
% =========================================================================
 
fig   = figure('Color','w','Position',[100 100 1280 1280]);
tl    = tiledlayout(fig, 3, 2);
 
axTop = nexttile(tl, [1 2]);   % Top panel  — full-width scene overview
axA   = nexttile(tl);          % Tile A     — corner 1 (FL by default ordering)
axB   = nexttile(tl);          % Tile B     — corner 3 (FR)
axC   = nexttile(tl);          % Tile C     — corner 2 (RR)
axD   = nexttile(tl);          % Tile D     — corner 4 (RL)
% =========================================================================
%% ── Section 3: Scenario Loading and Parameter Table Construction ─────────
% =========================================================================
 
% Load vehicle-location and pre-computed interference-path tables.
Scenario    = loadScenario(cfg.inputMat);
 
% Extract the unique set of vehicle IDs present in the scenario.
Vehicle_IDs = unique(double(Scenario.vLoc_table.Vid));
 
% Build the full radar-parameter table for all vehicles.
% Each row corresponds to one (vehicle, corner) sensor instance.
Parameters  = ARIFOMA_build_parameter_table(Vehicle_IDs, Scenario, cfg);
 
% =========================================================================
%% ── Section 4: Victim Vehicle Selection ─────────────────────────────────
% =========================================================================
 
if cfg.worstVehicle
    % Automated worst-case search: find the vehicle that experiences the
    % highest aggregate interference across all corners and snapshots.
    worst = ARIFOMA_find_worst_interference_vehicle( ...
                Scenario.vLoc_table, Scenario.intf_table, ...
                Parameters, cfg, opt);
    VehID = worst.vehIdWorst;
 
elseif cfg.specificVehicle
    % Use the fixed victim ID supplied in the configuration block.
    VehID = cfg.VehID;
 
else
    % Random selection: draw uniformly from sensor-equipped vehicles.
    validVeh = unique(Parameters.VehID(Parameters.HasSensor));
    VehID    = validVeh(randi(numel(validVeh)));
end
 
fprintf('Chosen vehicle: ID=%g\n', VehID);
 
% ── Apply victim-specific parameter overrides ─────────────────────────────
%   When the parameter mode supports it, the victim's waveform and scheme
%   are forced to the values in cfg.adaptWaveformVictim / cfg.adaptSchemeVictim.
isVict = (Parameters.VehID == VehID) & (Parameters.HasSensor);
if cfg.paramMode ~= "COMPASS"
    Parameters = ARIFOMA_apply_victim_overrides(Parameters, VehID, cfg);
end
 
% ── Victim-wide RF band summary ───────────────────────────────────────────
%   Used by downstream plotting routines and spectrum-allocation checks.
VictimVehicle = (Parameters.VehID == VehID) & Parameters.HasSensor;
bandMin_V     = min(Parameters.bandMin(VictimVehicle));
bandMax_V     = max(Parameters.bandMax(VictimVehicle));
 
% ── Geometry and visualisation options for frame extraction ───────────────
cfg.VehicleId       = VehID;
cfg.Role            = "victim";
cfg.Corner          = "all";
cfg.Tol             = 0.02;    %      Time-matching tolerance for frame extraction [s]
cfg.MaxPaths        = 300;     %      Maximum number of interference paths rendered
cfg.StretchY        = 0.4;     %      Vertical exaggeration factor for lane display
cfg.ShowAllVehicles = true;    %      Render all vehicles (not just victim + interferers)
cfg.WinX            = 150;     % [m]  Half-width of the scene window centred on victim
cfg.PlotGeometry    = true;    %      Draw path geometry overlays
cfg.MaxPathLen      = 150;     % [m]  Clip interference paths beyond this distance
cfg.FrontOnly       = false;   %      Include rear-facing corners in the view
 
% =========================================================================
%% ── Section 5: Radar Context Initialisation ──────────────────────────────
% =========================================================================
%   `radarCtx` is the central accumulator struct.  It is indexed by corner
%   name (e.g., radarCtx.FL) and holds per-CPI logs, scheme memory, and
%   frequency-history entries for every sensor encountered during the run.
% =========================================================================
 
radarCtx           = struct();
radarCtx.VehicleID = VehID;
radarCtx           = ARIFOMA_init_radarcorner_contexts(radarCtx, cfg.corners);
 
% Pre-populate the last-known carrier frequency (F0) for every sensor so
% that the first CPI step can reference a valid starting frequency without
% requiring special-case logic.
radarCtx = prepopulate_lastF0BySensor(radarCtx, cfg.corners, Parameters);
 
% =========================================================================
%% ── Section 6: Main Acquisition Time Loop ───────────────────────────────
% =========================================================================
%   Outer loop — iterates over acquisition time snapshots.
%   Inner loop — processes cfg.numberofCPIsperAcquisition CPI micro-steps
%                per snapshot for each of the four vehicle corners.
% =========================================================================
 
for it = 1 : numel(cfg.simTime)
 
    acquisition_t = cfg.simTime(it);   % Current acquisition time [s]
 
    % ── 6.1  Scene Snapshot Extraction and Top-Panel Rendering ───────────
    cla(axTop);
 
    % Extract a spatial snapshot of the scene at the current time.
    frame = ARIFOMA_frame_extract_mc( ...
                Scenario.vLoc_table, Scenario.intf_table, ...
                Parameters, acquisition_t, cfg);
 
    % Draw the global scene overview in the top panel (if enabled).
    if cfg.PlotGeometry
        ARIFOMA_frame_draw_mc(axTop, frame, cfg);
    end
 
    % ── 6.2  Build Per-Corner Link Tables and Scene Views ─────────────────
    %   `ARIFOMA_prepare_scene_views` constructs the per-corner
    %   victim / interferer views and the supporting link / path structures
    %   consumed by the CPI processing loop below.
    tolLinks = 0.002;   % [s] Tolerance for link-table time matching
    [cornerData, pathByCorner, links_sym, radarCtx] = ...
        ARIFOMA_prepare_scene_views( ...
            frame, Scenario.vLoc_table, Parameters, ...
            cfg.corners, radarCtx, tolLinks);
 
    % ── 6.3  CPI Micro-Step Loop ──────────────────────────────────────────
    %   Each CPI micro-step updates all four corner contexts independently.
    %   Corner processing order matches the axis assignment in Section 2.
    for singleCPI = 1 : cfg.numberofCPIsperAcquisition
 
        % Corner 1 — axA tile  (index 1 → FL)
        radarCtx = runCPI_perCorner(axA, 1, cornerData, radarCtx, ...
                                    pathByCorner, cfg, singleCPI, it, acquisition_t);
 
        % Corner 2 — axC tile  (index 2 → RR)
        radarCtx = runCPI_perCorner(axC, 2, cornerData, radarCtx, ...
                                    pathByCorner, cfg, singleCPI, it, acquisition_t);
 
        % Corner 3 — axD tile  (index 4 → RL)
        radarCtx = runCPI_perCorner(axD, 4, cornerData, radarCtx, ...
                                    pathByCorner, cfg, singleCPI, it, acquisition_t);
 
        % Corner 4 — axB tile  (index 3 → FR)
        radarCtx = runCPI_perCorner(axB, 3, cornerData, radarCtx, ...
                                    pathByCorner, cfg, singleCPI, it, acquisition_t);
 
    end  % singleCPI
end  % it
 
% Attach the final victim parameter snapshot to the context for downstream
% post-processing scripts that need corner geometry without reloading the
% full scenario.
radarCtx.victimData = cornerData.victim;
 
% =========================================================================
%% ── Section 7: Save Results ──────────────────────────────────────────────
% =========================================================================
 
outFile = sprintf('radarCtx%s.mat', cfg.sceneTag);
save(outFile, "radarCtx", "-v7.3");
fprintf('Results saved to: %s\n', outFile);
% =========================================================================
%% ── Local Helper Functions ───────────────────────────────────────────────
% =========================================================================
 
% -------------------------------------------------------------------------
function S = loadScenario(matPath)
% LOADSCENARIO  Load the required scenario variables from a MAT-file.
%
%   S = LOADSCENARIO(matPath) reads the MAT-file located at matPath and
%   returns a struct S containing the two tables that ARIFOMA requires:
%
%     S.vLoc_table — vehicle location table (one row per vehicle × time)
%     S.intf_table — pre-computed interference-path table
%
%   The function calls error() immediately if either variable is absent,
%   preventing silent failures later in the simulation pipeline.
%
%   INPUT
%     matPath  (char | string) — absolute or relative path to the MAT-file.
%
%   OUTPUT
%     S  (struct) — workspace struct loaded from matPath, with the two
%                   required fields guaranteed to be present.
 
    S = load(matPath);
 
    % Validate vLoc_table presence.
    try
        vLoc_table = S.vLoc_table;                           %#ok<NASGU>
    catch
        error('ARIFOMA:loadScenario:missingVariable', ...
              'Input MAT-file must contain ''vLoc_table''.');
    end
 
    % Validate intf_table presence.
    try
        intf_table = S.intf_table;                           %#ok<NASGU>
    catch
        error('ARIFOMA:loadScenario:missingVariable', ...
              'Input MAT-file must contain ''intf_table''.');
    end
 
    S.vLoc_table = S.vLoc_table;
    S.intf_table = S.intf_table;
end
 
% -------------------------------------------------------------------------
function radarCtx = runCPI_perCorner( ...
        ax, cIndex, cornerData, radarCtx, pathByCorner, cfg, ...
        numCPI, tIndex, t0_cur)
% RUNCPI_PERCORNER  Execute one CPI micro-step for a single victim corner.
%
%   radarCtx = RUNCPI_PERCORNER(ax, cIndex, cornerData, radarCtx,
%   pathByCorner, cfg, numCPI, tIndex, t0_cur) performs the full
%   per-corner CPI processing pipeline for corner index cIndex and stores
%   the updated context back into radarCtx.
%
%   PIPELINE STAGES
%     1. Select the corner-specific victim / interferer parameter bundle.
%     2. Prepare the target axis and return early if no sensor is present.
%     3. Apply scheme logic (frequency hopping, BATS, etc.) for all sensors.
%     4. Build the engine input bundle consumed by the strategy dispatcher.
%     5. Dispatch to the selected waveform / mitigation strategy.
%     6. Post-process output metrics and persist the updated corner context.
%
%   INPUTS
%     ax          (Axes)    — target tile axes for visualisation.
%     cIndex      (int)     — corner index into cornerData (1–4).
%     cornerData  (struct[])— per-corner scene bundle from prepare_scene_views.
%     radarCtx    (struct)  — accumulated radar context struct.
%     pathByCorner(struct)  — per-corner interference path structures.
%     cfg         (struct)  — global configuration struct.
%     numCPI      (int)     — current CPI micro-step index (1-based).
%     tIndex      (int)     — current acquisition time index (1-based).
%     t0_cur      (double)  — current acquisition time [s].
%
%   OUTPUT
%     radarCtx  — updated struct with the processed corner context stored
%                 under radarCtx.<cornerName>.
 
    numberofCPIsperAcquisition = cfg.numberofCPIsperAcquisition;
    targetSpan                 = cfg.CPI_duration_s;
 
    % ── Stage 1: Select corner-specific scene objects ─────────────────────
    corner = get_corner_inputs(cIndex, cornerData, radarCtx);
    cn     = corner.cn;
    ctx    = corner.ctx;
 
    % ── Stage 2: Axis preparation + early exit if no sensor present ───────
    prep_axes(ax);
    if isempty(corner.V)
        title(ax, sprintf('%s: no sensor', cn));
        radarCtx.(cn) = ctx;
        return;
    end
 
    % ── Stage 3: Apply per-sensor scheme logic ────────────────────────────
    %   Updates carrier frequency, hop state, and sensing memory for both
    %   the victim and all active interferers in this corner view.
    [corner, radarCtx] = ARIFOMA_apply_all_schemes( ...
        corner, radarCtx, numCPI, numberofCPIsperAcquisition, tIndex);
 
    % ── Stage 4: Build engine input bundle ───────────────────────────────
    engine = ARIFOMA_build_engine_inputs( ...
        ax, corner, pathByCorner, cfg, numCPI, ...
        numberofCPIsperAcquisition, tIndex, t0_cur);
 
    % ── Stage 5: Dispatch to selected strategy ────────────────────────────
    out = ARIFOMA_run_strategy(engine);
 
    % ── Stage 6: Post-process metrics and persist context ─────────────────
    ctx = ARIFOMA_postprocess_and_update_ctx( ...
        corner.ctx, out, corner, cfg, numCPI, ...
        numberofCPIsperAcquisition, tIndex, t0_cur);
 
    radarCtx.(cn) = ctx;
end
 
% -------------------------------------------------------------------------
function corner = get_corner_inputs(cIndex, cornerData, radarCtx)
% GET_CORNER_INPUTS  Extract the corner-specific victim / interferer bundle.
%
%   corner = GET_CORNER_INPUTS(cIndex, cornerData, radarCtx) builds a
%   lightweight struct that aggregates all objects needed for one corner's
%   CPI processing step.
%
%   INPUTS
%     cIndex     (int)     — 1-based corner index into cornerData.
%     cornerData (struct[])— array of per-corner scene bundles.
%     radarCtx   (struct)  — accumulated radar context.
%
%   OUTPUT
%     corner (struct) with fields:
%       .cn           — corner name string (e.g., 'FL')
%       .V            — victim sensor parameter row
%       .I            — interferer parameter rows
%       .perspectives — per-sensor perspective structs
%       .pairsBySensor— per-sensor interferer pairing map
%       .ctx          — current corner context from radarCtx
 
    corner.cn            = char(cornerData(cIndex).corner);
    corner.V             = cornerData(cIndex).victim;
    corner.I             = cornerData(cIndex).interferers;
    corner.perspectives  = cornerData(cIndex).perspective;
    corner.pairsBySensor = cornerData(cIndex).pairsBySensor;
    corner.ctx           = radarCtx.(corner.cn);
end