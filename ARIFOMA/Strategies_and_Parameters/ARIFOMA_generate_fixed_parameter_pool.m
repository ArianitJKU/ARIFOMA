function Param = ARIFOMA_generate_fixed_parameter_pool(vehicleIds, cfg)
% ARIFOMA_GENERATE_FIXED_POOL_PARAMETERS
%
% Generates a fixed pool of cfg.numFixedWaveforms unique FMCW radar
% parameter configurations and assigns exactly one configuration per
% vehicle, applied uniformly across all corners of that vehicle.
%
% This enforces a constrained deployment where the total number of
% distinct waveform types in the scene is limited to numFixedWaveforms,
% regardless of the number of vehicles. Each vehicle draws once from the
% pool at random; all corners of that vehicle inherit the same waveform.
%
% INPUTS
%   vehicleIds  (double array) — unique vehicle IDs in the scenario
%   cfg         (struct)       — global configuration struct; must contain:
%                                  cfg.numFixedWaveforms  (int)
%                                  cfg.corners            (string array)
%                                  cfg.CPI_duration_s     (double)
%                                  cfg.pHasSensor         (double)
%                                  cfg.UseFixedSeed       (logical)
%                                  cfg.Seed               (double)
%                                  cfg.bandByScenario     (struct)
%                                  cfg.scenarioDist       (struct)
%                                  cfg.SchemeList/Weights (string/double)
%                                  cfg.WaveformList/Weights
%                                  cfg.StepCountMin/Max
%                                  cfg.StepHzMin/Max
%                                  cfg.FlybackDefaultMin/Max_s
%                                  cfg.FlybackSteppedMin/Max_s
%                                  cfg.ReturnFlyMin/Max_s
%                                  cfg.SensorIdMap        (table)
%                                  cfg.KeepIdForAbsent    (logical)
%
% OUTPUT
%   Param  (table) — parameter table with one row per (vehicle, corner)
%                    pair, matching the schema of
%                    ARIFOMA_generate_corner_parameters.

%% ── RNG Initialisation ───────────────────────────────────────────────────
if cfg.UseFixedSeed
    rng(cfg.Seed);
end

veh     = vehicleIds(:);
corners = string(cfg.corners(:));
nV      = numel(veh);
nC      = numel(corners);
nPool   = cfg.numFixedWaveforms;
PRI_s   = cfg.CPI_duration_s;

%% ── Step 1: Generate the Fixed Waveform Pool ─────────────────────────────
%   Create exactly nPool unique parameter configurations.
%   These are the only waveform types that will exist in the scene.

pool = generate_pool(nPool, cfg, PRI_s);

fprintf('[ARIFOMA] Fixed pool: %d unique waveform configurations generated.\n', nPool);

%% ── Step 2: Assign One Pool Entry per Vehicle ────────────────────────────
%   Each vehicle draws one index from the pool at random.
%   All corners of that vehicle inherit the same configuration.

vehiclePoolIndex = randi(nPool, nV, 1);   % one pool index per vehicle

%% ── Step 3: Build the Parameter Table ───────────────────────────────────
nRows = nV * nC;

S = struct( ...
    'VehID',[], 'Corner',[], 'SensorID',[], 'HasSensor',[], ...
    'Scenario',[], 'Scheme',[], 'Waveform',[], 'WaveformClean',[], 'Dir',[], ...
    'B',[], 'Tch',[], 'Nch',[], 'ChirpsPerStep',[], 'FreqIncrHz',[], ...
    'F0',[], 'Fc',[], ...
    'T_flyback_default',[], 'T_flyback_stepped',[], 'T_returnfly',[], ...
    'bandMin',[], 'bandMax',[], ...
    'Nseg',[], 'PermPattern',[]);

S = repmat(S, nRows, 1);
row = 0;

for iV = 1:nV

    VehID         = veh(iV);
    vehHasSensors = (rand <= cfg.pHasSensor);
    p             = pool(vehiclePoolIndex(iV));   % this vehicle's waveform config

    for iC = 1:nC
        row    = row + 1;
        Corner = corners(iC);

        % ── Sensor ID ────────────────────────────────────────────────────
        SensorID = VehID * 100 + corner_index(Corner);

        if ~isempty(cfg.SensorIdMap)
            hit = cfg.SensorIdMap( ...
                cfg.SensorIdMap.VehID   == VehID & ...
                cfg.SensorIdMap.Corner  == Corner, :);
            if ~isempty(hit)
                SensorID = double(hit.SensorID(1));
            end
        end

        if ~vehHasSensors && ~cfg.KeepIdForAbsent
            SensorID = NaN;
        end

        % ── Segmentation (sensor-specific, deterministic) ─────────────
        if isfinite(SensorID)
            Nseg = 4 + mod(SensorID, 5);
            rs   = RandStream('mt19937ar', 'Seed', uint32(2654435761*SensorID + iC));
            permPattern = randperm(rs, Nseg);
        else
            Nseg        = 4;
            permPattern = 1:4;
        end

        % ── Assign pool configuration or empty defaults ───────────────
        if vehHasSensors
            S(row) = struct( ...
                'VehID',          VehID, ...
                'Corner',         Corner, ...
                'SensorID',       SensorID, ...
                'HasSensor',      vehHasSensors, ...
                'Scenario',       p.Scenario, ...
                'Scheme',         p.Scheme, ...
                'Waveform',       p.Waveform, ...
                'WaveformClean',  p.WaveformClean, ...
                'Dir',            p.Dir, ...
                'B',              p.B, ...
                'Tch',            p.Tch, ...
                'Nch',            p.Nch, ...
                'ChirpsPerStep',  p.ChirpsPerStep, ...
                'FreqIncrHz',     p.FreqIncrHz, ...
                'F0',             p.F0, ...
                'Fc',             p.Fc, ...
                'T_flyback_default', p.T_flyback_default, ...
                'T_flyback_stepped', p.T_flyback_stepped, ...
                'T_returnfly',    p.T_returnfly, ...
                'bandMin',        p.bandMin, ...
                'bandMax',        p.bandMax, ...
                'Nseg',           Nseg, ...
                'PermPattern',    permPattern);
        else
            % Vehicle has no sensor — fill with NaN sentinels
            S(row) = struct( ...
                'VehID',          VehID, ...
                'Corner',         Corner, ...
                'SensorID',       SensorID, ...
                'HasSensor',      false, ...
                'Scenario',       "", ...
                'Scheme',         "", ...
                'Waveform',       "", ...
                'WaveformClean',  "", ...
                'Dir',            NaN, ...
                'B',              NaN, ...
                'Tch',            NaN, ...
                'Nch',            NaN, ...
                'ChirpsPerStep',  NaN, ...
                'FreqIncrHz',     NaN, ...
                'F0',             NaN, ...
                'Fc',             NaN, ...
                'T_flyback_default', NaN, ...
                'T_flyback_stepped', NaN, ...
                'T_returnfly',    NaN, ...
                'bandMin',        NaN, ...
                'bandMax',        NaN, ...
                'Nseg',           Nseg, ...
                'PermPattern',    permPattern);
        end
    end
end

Param = struct2table(S);

%% ── Verification Report ──────────────────────────────────────────────────
%   Confirm that the number of unique F0 values does not exceed nPool.
uniqueF0 = numel(unique(Param.F0(Param.HasSensor)));
fprintf('[ARIFOMA] Unique F0 values in scene: %d (pool size: %d)\n', ...
        uniqueF0, nPool);
if uniqueF0 > nPool
    warning('ARIFOMA:poolOverflow', ...
        ['More unique F0 values than pool size detected. ' ...
         'Check pool generation and assignment logic.']);
end

end

%% =========================================================================
%% Local Functions
%% =========================================================================

function pool = generate_pool(nPool, cfg, PRI_s)
% GENERATE_POOL  Generate nPool unique waveform parameter configurations.
%
%   Each entry is drawn randomly from the configured scenario, scheme,
%   waveform, and RF band distributions. The pool is generated once and
%   then reused for all vehicle assignments.

    pool = struct( ...
        'Scenario',[], 'Scheme',[], 'Waveform',[], 'WaveformClean',[], 'Dir',[], ...
        'B',[], 'Tch',[], 'Nch',[], 'ChirpsPerStep',[], 'FreqIncrHz',[], ...
        'F0',[], 'Fc',[], ...
        'T_flyback_default',[], 'T_flyback_stepped',[], 'T_returnfly',[], ...
        'bandMin',[], 'bandMax',[]);

    pool = repmat(pool, nPool, 1);

    for p = 1:nPool

        % Draw scenario from a random corner reference distribution
        cornerRef = cfg.corners(randi(numel(cfg.corners)));
        Scenario  = pick_weighted( ...
            cfg.scenarioDist.(cornerRef).list, ...
            cfg.scenarioDist.(cornerRef).weights);

        sc      = cfg.bandByScenario.(Scenario);
        bandMin = sc.bandMin;
        bandMax = sc.bandMax;
        Brange  = sc.Brange;
        Trange  = sc.Trange;

        Scheme   = pick_weighted(cfg.SchemeList,   cfg.SchemeWeights);
        Waveform = pick_weighted(cfg.WaveformList, cfg.WaveformWeights);
        [WaveformClean, Dir] = normalize_Waveform(Waveform);

        B   = Brange(1) + rand*(Brange(2) - Brange(1));
        Tch = Trange(1) + rand*(Trange(2) - Trange(1));

        [lo, hi] = order_pair(cfg.FlybackDefaultMin_s, cfg.FlybackDefaultMax_s);
        Tfly_def  = pick_real(lo, hi);

        [lo, hi] = order_pair(cfg.FlybackSteppedMin_s, cfg.FlybackSteppedMax_s);
        Tfly_step = pick_real(lo, hi);

        [lo, hi] = order_pair(cfg.ReturnFlyMin_s, cfg.ReturnFlyMax_s);
        Tret      = pick_real(lo, hi);

        fly_s = pick_flyback(WaveformClean, Tfly_def, Tfly_step);
        Nch   = max(1, floor(PRI_s / max(eps, Tch + fly_s)));

        ChirpsPerStep = NaN;
        FreqIncrHz    = 0;
        nStepForF0    = 1;

        if WaveformClean == "STEPPED FMCW"
            [lo, hi] = order_pair(cfg.StepCountMin, cfg.StepCountMax);
            Nsteps = pick_int(max(1,floor(lo)), max(1,floor(hi)));
            Nsteps = min(Nsteps, Nch);

            [lo, hi] = order_pair(cfg.StepHzMin, cfg.StepHzMax);
            stepHz = pick_real(lo, hi);

            [B, stepHz]   = fit_span_to_band(bandMin, bandMax, B, Nsteps, stepHz);
            ChirpsPerStep = Nsteps;
            FreqIncrHz    = stepHz;
            nStepForF0    = Nsteps;
        end

        [F0, Fc] = choose_f0_within_band( ...
            bandMin, bandMax, B, nStepForF0, Dir, FreqIncrHz);

        pool(p) = struct( ...
            'Scenario',       Scenario, ...
            'Scheme',         Scheme, ...
            'Waveform',       Waveform, ...
            'WaveformClean',  WaveformClean, ...
            'Dir',            Dir, ...
            'B',              B, ...
            'Tch',            Tch, ...
            'Nch',            Nch, ...
            'ChirpsPerStep',  ChirpsPerStep, ...
            'FreqIncrHz',     FreqIncrHz, ...
            'F0',             F0, ...
            'Fc',             Fc, ...
            'T_flyback_default', Tfly_def, ...
            'T_flyback_stepped', Tfly_step, ...
            'T_returnfly',    Tret, ...
            'bandMin',        bandMin, ...
            'bandMax',        bandMax);
    end
end

% ── Inline helpers (mirrored from ARIFOMA_generate_corner_parameters) ─────

function fly_s = pick_flyback(WaveformClean, Tfly_def, Tfly_step)
    if WaveformClean == "STEPPED FMCW"
        fly_s = Tfly_step;
    else
        fly_s = Tfly_def;
    end
end

function [B_adj, stepHz_adj] = fit_span_to_band(bandMin, bandMax, B, Nch, stepHz)
    bandMin = double(bandMin); bandMax = double(bandMax);
    B = double(B); Nch = double(Nch); stepHz = double(stepHz);
    bwAvail = bandMax - bandMin;
    if Nch <= 1
        B_adj = min(B, bwAvail); stepHz_adj = 0; return;
    end
    if B > bwAvail, B = bwAvail; end
    maxStepHz  = max(0, (bwAvail - B) / (Nch - 1));
    stepHz_adj = min(max(0, stepHz), maxStepHz);
    span = B + (Nch-1)*stepHz_adj;
    if span > bwAvail
        B = max(0, bwAvail - (Nch-1)*stepHz_adj);
    end
    B_adj = B;
end

function [wfClean, dir] = normalize_Waveform(label, dirExisting)
    L = upper(string(label));
    switch L
        case "POSITIVE-SLOPE",         wfClean = "FMCW";         dir = +1;
        case "NEGATIVE-SLOPE",         wfClean = "FMCW";         dir = -1;
        case "POSITIVE-SLOPE STEPPED", wfClean = "STEPPED FMCW"; dir = +1;
        case "NEGATIVE-SLOPE STEPPED", wfClean = "STEPPED FMCW"; dir = -1;
        otherwise,                     wfClean = "FMCW";         dir = +1;
    end
    if nargin >= 2 && isfinite(dirExisting) && dirExisting ~= 0
        dir = sign(dirExisting);
    end
end

function x = pick_weighted(values, weights)
    values = string(values(:));
    w = weights(:) / sum(weights(:));
    x = values(find(rand <= cumsum(w), 1, 'first'));
end

function [F0, Fc] = choose_f0_within_band(bandMin, bandMax, B, Nch, dir, stepHz)
    if ~isfinite(stepHz) || stepHz <= 0, stepHz = 0; end
    span = max(B + (Nch-1)*stepHz, B);
    if span > (bandMax - bandMin)
        B    = min(B, (bandMax-bandMin) - max(0,(Nch-1)*stepHz));
        span = B + (Nch-1)*stepHz;
    end
    base = bandMin + rand * max(0, (bandMax - bandMin) - span);
    F0   = base + (dir < 0)*B;
    Fc   = base + B/2;
end

function [lo, hi] = order_pair(a, b)
    lo = min(double(a), double(b));
    hi = max(double(a), double(b));
end

function x = pick_real(lo, hi)
    x = lo + (hi > lo) * rand * (hi - lo);
end

function k = pick_int(lo, hi)
    if lo == hi, k = lo; else, k = randi([lo hi]); end
end

function idx = corner_index(Corner)
    map = struct('FL',1,'FR',2,'RL',3,'RR',4,'F',5);
    if isfield(map, char(Corner))
        idx = map.(char(Corner));
    else
        idx = 0;
    end
end