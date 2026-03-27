%% main_sim_FF.m — ARIFOMA (Front-Facing Sensor)
% =========================================================================
%
% Automotive Radar Interference Figures of Merit Analysis Framework
%
% PURPOSE
% -------
%   Simulates FMCW radar-level interference on a multi-lane highway
%   scene and evaluates interference mitigation strategies for a single
%   front-facing (FF) automotive radar sensor.
%
% SIMULATION PIPELINE
% -------------------
%   1. Load a pre-computed traffic / propagation scenario from disk.
%   2. Generate radar parameters for all sensor-equipped vehicles.
%   3. Select a victim vehicle (worst-case, specific ID, or random).
%   4. Iterate over acquisition time snapshots (cfg.simTime).
%      a. Extract and draw the per-snapshot scene geometry.
%      b. Build the front-facing (FF) victim / interferer link table.
%      c. For each CPI micro-step, run the selected waveform / mitigation
%         strategy and accumulate the resulting radar context.
%   5. Save the final `radarCtx` struct for post-processing and analysis.
%
% OUTPUT
% ------
%   radarCtx<cfg.sceneTag>.mat — FF radar context struct containing
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
%   • This variant uses a single front-facing corner labelled "FF".
%   • All corner-indexed logic from the four-corner version is reduced
%     to operate on the single "FF" entry only.
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
cfg.FrontSensor = true;
% ── 1.1  Simulation Time Axis ─────────────────────────────────────────────
cfg.TA      = 0.1;
cfg.TS      = 30;
cfg.simTime = 0:cfg.TA:cfg.TS;          % [s]  Acquisition times (0 – 30 s, 0.1 s steps)

% ── 1.2  Scenario / Output Naming ────────────────────────────────────────
cfg.sceneTag    = '_FF_SRR_25MHz_density85';
cfg.bats_eps    = 0.07;                 %      BATS algorithm convergence threshold
cfg.outVideoFile = sprintf('scene%s.mp4', cfg.sceneTag);

% ── 1.3  Playback / Animation Settings ───────────────────────────────────
cfg.fps            = 24;                % [fps] Output video frame rate
cfg.holdTopFrames  = 1;                 % [frames] Hold count for top-panel frame
cfg.holdTileFrames = 2;                 % [frames] Hold count for FF-tile frame

% ── 1.4  Chirp Flyback / Return-fly Timing Limits ────────────────────────
cfg.FlybackDefaultMin_s = 3e-6;
cfg.FlybackDefaultMax_s = 10e-6;
cfg.FlybackSteppedMin_s = 2e-6;
cfg.FlybackSteppedMax_s = 4e-6;
cfg.ReturnFlyMin_s      = 5e-6;
cfg.ReturnFlyMax_s      = 7e-6;

% ── 1.5  Per-Class RF Band and Waveform Parameter Ranges ─────────────────
cfg.bandByScenario = struct();

cfg.bandByScenario.MRR = struct( ...
    'bandMin',           77e9, ...
    'bandMax',           81e9, ...
    'Brange',     [400e6 800e6], ...
    'Trange',    [20e-6  40e-6], ...
    'T_flyback_default', [3e-6 5e-6], ...
    'T_flyback_stepped', [1e-6 3e-6], ...
    'T_returnfly',       [4e-6 6e-6]);

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
cfg.paramMode        = "EXTENDED";
cfg.UseFixedSeed     = false;
cfg.numFixedWaveforms = 10;
cfg.Seed             = 7;
cfg.pHasSensor       = 1;

% ── 1.7  Victim Corner Scenario Assignment ────────────────────────────────
%   Single front-facing sensor — only FF corner is defined.
cfg.victimscenarioPerCorner = struct('FF', "MRR");

% ── 1.8  Radar-Class Distribution for the General Traffic Population ──────
cfg.scenarioDist.FF.list    = ["MRR","SRR","LRR"];
cfg.scenarioDist.FF.weights = [1 0 0];

% ── 1.9  Waveform and Scheme Pools for Traffic Population ─────────────────
cfg.SchemeList    = ["CONVENTIONAL","RANDOM_FH_PER_CPI","FAR", ...
                     "BLUEFMCW","BATSINSPIRED_FH_PER_CPI"];
cfg.SchemeWeights = [1 0 0 0 0];

cfg.WaveformList    = ["POSITIVE-SLOPE","POSITIVE-SLOPE STEPPED", ...
                       "NEGATIVE-SLOPE","NEGATIVE-SLOPE STEPPED"];
cfg.WaveformWeights = [1 0 0 0];

cfg.StepCountMin = 64;
cfg.StepCountMax = 128;
cfg.StepHzMin    = 10e6;
cfg.StepHzMax    = 20e6;

% ── 1.10  ADC and Detector Settings ──────────────────────────────────────
cfg.BADC_Hz       = 25e6;
cfg.randstartTime = 40e-6;

% ── 1.11  COMPASS Mode: Lane and Band Settings ────────────────────────────
cfg.numberoflanes = 6;
cfg.lane_edges    = 4 + (0 : cfg.numberoflanes-1) * 4;
cfg.compassShare  = 1;

cfg.CPI_duration_s = 10e-3;

% ── 1.12  Plotting and Debug Flags ───────────────────────────────────────
cfg.plotScenario                = false;
cfg.plotScenarioEnd             = true;
cfg.boolIgnoreFlybackInterferer = true;

% ── 1.13  COMPASS Carrier-Frequency Band Boundaries ──────────────────────
cfg.compassFcBands.NW = [77e9 78e9];
cfg.compassFcBands.NE = [78e9 79e9];
cfg.compassFcBands.SE = [79e9 80e9];
cfg.compassFcBands.SW = [80e9 81e9];

% ── 1.14  Victim Vehicle Selection ───────────────────────────────────────
cfg.worstVehicle    = false;
cfg.specificVehicle = true;
cfg.VehID           = 243;

% ── 1.15  Victim Waveform / Scheme Override ───────────────────────────────
cfg.adaptWaveformVictim = "POSITIVE-SLOPE";
cfg.adaptSchemeVictim   = "BATSINSPIRED_FH_PER_CPI";

% ── 1.16  CPI Micro-Steps per Acquisition Snapshot ───────────────────────
cfg.numberofCPIsperAcquisition = 2;

% ── 1.17  Corner Ordering ─────────────────────────────────────────────────
%   Single front-facing sensor — only "FF" corner.
cfg.corners     = ["FF"];
cfg.SensorIdMap = table();
cfg.KeepIdForAbsent = true;

% ── 1.18  Input Scenario File ─────────────────────────────────────────────
cfg.inputMat = 'C:\Users\PC\Desktop\ARIFOMA\SUMO_WiLabVISIM_Data\Highway_Density85_3Lanes.mat';

% =========================================================================
%% ── Section 2: Figure and Tiled-Layout Initialisation ───────────────────
% =========================================================================
%   A 2×1 tiled layout is created: top panel for the global scene overview,
%   bottom panel for the single FF sensor view.
% =========================================================================

fig   = figure('Color','w','Position',[100 100 1280 860]);
tl    = tiledlayout(fig, 2, 1);

axTop = nexttile(tl);   % Top panel — full-width scene overview
axFF  = nexttile(tl);   % Bottom panel — FF sensor view

% =========================================================================
%% ── Section 3: Scenario Loading and Parameter Table Construction ─────────
% =========================================================================

Scenario    = loadScenario(cfg.inputMat);
Vehicle_IDs = unique(double(Scenario.vLoc_table.Vid));
Parameters  = ARIFOMA_build_parameter_table(Vehicle_IDs, Scenario, cfg);

% =========================================================================
%% ── Section 4: Victim Vehicle Selection ─────────────────────────────────
% =========================================================================

if cfg.worstVehicle
    worst = ARIFOMA_find_worst_interference_vehicle( ...
                Scenario.vLoc_table, Scenario.intf_table, ...
                Parameters, cfg, opt);
    VehID = worst.vehIdWorst;

elseif cfg.specificVehicle
    VehID = cfg.VehID;

else
    validVeh = unique(Parameters.VehID(Parameters.HasSensor));
    VehID    = validVeh(randi(numel(validVeh)));
end

fprintf('Chosen vehicle: ID=%g\n', VehID);

% ── Apply victim-specific parameter overrides ─────────────────────────────
if cfg.paramMode ~= "COMPASS"
    Parameters = ARIFOMA_apply_victim_overrides(Parameters, VehID, cfg);
end

% ── Victim-wide RF band summary ───────────────────────────────────────────
VictimVehicle = (Parameters.VehID == VehID) & Parameters.HasSensor;
bandMin_V     = min(Parameters.bandMin(VictimVehicle));
bandMax_V     = max(Parameters.bandMax(VictimVehicle));

% ── Geometry and visualisation options ───────────────────────────────────
cfg.VehicleId       = VehID;
cfg.Role            = "victim";
cfg.Corner          = "FF";
cfg.Tol             = 0.02;
cfg.MaxPaths        = 300;
cfg.StretchY        = 0.4;
cfg.ShowAllVehicles = true;
cfg.WinX            = 150;
cfg.PlotGeometry    = true;
cfg.MaxPathLen      = 150;
cfg.FrontOnly       = true;    % Front-facing only — rear corners excluded

% =========================================================================
%% ── Section 5: Radar Context Initialisation ──────────────────────────────
% =========================================================================

radarCtx           = struct();
radarCtx.VehicleID = VehID;
radarCtx           = ARIFOMA_init_radarcorner_contexts(radarCtx, cfg.corners);
radarCtx           = prepopulate_lastF0BySensor(radarCtx, cfg.corners, Parameters);

% =========================================================================
%% ── Section 6: Main Acquisition Time Loop ───────────────────────────────
% =========================================================================

for it = 1 : numel(cfg.simTime)

    acquisition_t = cfg.simTime(it);

    % ── 6.1  Scene Snapshot Extraction and Top-Panel Rendering ───────────
    cla(axTop);

    frame = ARIFOMA_frame_extract_mc( ...
                Scenario.vLoc_table, Scenario.intf_table, ...
                Parameters, acquisition_t, cfg);

    if cfg.PlotGeometry
        ARIFOMA_frame_draw_mc(axTop, frame, cfg);
    end

    % ── 6.2  Build FF Corner Link Table and Scene View ────────────────────
    tolLinks = 0.002;
    [cornerData, pathByCorner, links_sym, radarCtx] = ...
        ARIFOMA_prepare_scene_views( ...
            frame, Scenario.vLoc_table, Parameters, ...
            cfg.corners, radarCtx, tolLinks);

    % ── 6.3  CPI Micro-Step Loop ──────────────────────────────────────────
    for singleCPI = 1 : cfg.numberofCPIsperAcquisition

        % Single FF corner
        radarCtx = runCPI_perCorner(axFF, 1, cornerData, radarCtx, ...
                                    pathByCorner, cfg, singleCPI, it, acquisition_t);

    end  % singleCPI
end  % it

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
%   S = LOADSCENARIO(matPath) reads the MAT-file at matPath and returns a
%   struct S with the two tables ARIFOMA requires:
%
%     S.vLoc_table — vehicle location table (one row per vehicle × time)
%     S.intf_table — pre-computed interference-path table

    S = load(matPath);

    try
        vLoc_table = S.vLoc_table;                           %#ok<NASGU>
    catch
        error('ARIFOMA:loadScenario:missingVariable', ...
              'Input MAT-file must contain ''vLoc_table''.');
    end

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
% RUNCPI_PERCORNER  Execute one CPI micro-step for the FF sensor corner.
%
%   Identical pipeline to the four-corner version but operates solely on
%   the single FF corner entry in cornerData.

    numberofCPIsperAcquisition = cfg.numberofCPIsperAcquisition;

    % ── Stage 1: Select FF corner scene objects ───────────────────────────
    corner = get_corner_inputs(cIndex, cornerData, radarCtx);
    cn     = corner.cn;   % Will be 'FF'
    ctx    = corner.ctx;

    % ── Stage 2: Axis preparation + early exit if no sensor present ───────
    prep_axes(ax);
    if isempty(corner.V)
        title(ax, sprintf('%s: no sensor', cn));
        radarCtx.(cn) = ctx;
        return;
    end

    % ── Stage 3: Apply per-sensor scheme logic ────────────────────────────
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
% GET_CORNER_INPUTS  Extract the FF corner victim / interferer bundle.

    corner.cn            = char(cornerData(cIndex).corner);  % 'FF'
    corner.V             = cornerData(cIndex).victim;
    corner.I             = cornerData(cIndex).interferers;
    corner.perspectives  = cornerData(cIndex).perspective;
    corner.pairsBySensor = cornerData(cIndex).pairsBySensor;
    corner.ctx           = radarCtx.(corner.cn);
end