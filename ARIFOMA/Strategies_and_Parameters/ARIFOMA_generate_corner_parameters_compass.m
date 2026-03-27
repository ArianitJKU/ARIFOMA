function Param = ARIFOMA_generate_corner_parameters_compass(vehicleIds, cfg)
% ARIFOMA_GENERATE_CORNER_PARAMETERS_COMPASS
%   Hybrid COMPASS / random per-corner FMCW radar parameter generator.
%
% PURPOSE
% -------
%   Assigns FMCW waveform parameters to every (vehicle, corner) sensor
%   instance in the scene using a hybrid frequency-assignment strategy:
%
%     - The victim vehicle (cfg.VehID) always receives compass-based
%       carrier frequency assignment.
%     - A configurable fraction (cfg.compassShare) of the remaining
%       traffic population also receives compass assignment.
%     - All other vehicles receive fully random carrier frequency
%       assignment within the RF band (identical to EXTENDED mode).
%
%   In compass mode, each corner's carrier frequency is anchored to the
%   midpoint of a pre-defined 1 GHz sub-band determined by the corner's
%   effective facing direction (NW / NE / SE / SW), after applying a
%   180-degree inversion for vehicles in opposing traffic lanes.
%
% HYBRID MODE SUMMARY
% -------------------
%   Vehicle type          | Frequency assignment
%   ----------------------|---------------------------------------------
%   Victim (cfg.VehID)    | Always compass  (mandatory)
%   compassShare fraction | Compass         (random subset of population)
%   Remaining population  | Random within RF band (EXTENDED behaviour)
%
% INPUTS
%   vehicleIds  (double array) — column vector of unique vehicle IDs
%   cfg         (struct)       — global configuration struct; required fields:
%
%     RNG
%       cfg.UseFixedSeed        (logical) — seed the RNG before generation
%       cfg.Seed                (double)  — RNG seed value
%
%     Geometry
%       cfg.corners             (string array) — corner labels
%       cfg.SensorIdMap         (table)   — optional external sensor-ID lookup
%       cfg.KeepIdForAbsent     (logical) — retain IDs for sensor-less corners
%       cfg.pHasSensor          (double)  — probability a vehicle carries sensors
%       cfg.CornerInvertMap     (table)   — per-vehicle inversion flags
%                                           {VehID, Invert}; produced by
%                                           ARIFOMA_build_parameter_table
%
%     Timing
%       cfg.CPI_duration_s      (double)  — full PRI / CPI duration [s]
%       cfg.FlybackDefaultMin_s (double)  — default flyback lower bound [s]
%       cfg.FlybackDefaultMax_s (double)  — default flyback upper bound [s]
%       cfg.FlybackSteppedMin_s (double)  — stepped flyback lower bound [s]
%       cfg.FlybackSteppedMax_s (double)  — stepped flyback upper bound [s]
%       cfg.ReturnFlyMin_s      (double)  — return-fly lower bound [s]
%       cfg.ReturnFlyMax_s      (double)  — return-fly upper bound [s]
%
%     RF Band and Waveform
%       cfg.scenarioDist        (struct)  — per-corner radar-class distributions
%       cfg.bandByScenario      (struct)  — per-class RF and waveform ranges
%       cfg.SchemeList          (string array) — mitigation scheme pool
%       cfg.SchemeWeights       (double array) — scheme selection weights
%       cfg.WaveformList        (string array) — waveform type pool
%       cfg.WaveformWeights     (double array) — waveform selection weights
%
%     Stepped FMCW
%       cfg.StepCountMin        (double)  — minimum sub-sweep step count
%       cfg.StepCountMax        (double)  — maximum sub-sweep step count
%       cfg.StepHzMin           (double)  — minimum frequency step size [Hz]
%       cfg.StepHzMax           (double)  — maximum frequency step size [Hz]
%
%     COMPASS
%       cfg.VehID               (double)  — victim vehicle ID (always compass)
%       cfg.compassShare        (double)  — fraction of non-victim vehicles
%                                           assigned compass [0, 1]
%       cfg.compassFcBands      (struct)  — per-sector carrier-frequency bands
%                                           fields: NW, NE, SE, SW
%                                           each a [fLow fHigh] double pair [Hz]
%
% OUTPUT
%   Param  (table) — parameter table with one row per (vehicle, corner).
%
%   Additional columns vs. ARIFOMA_generate_corner_parameters:
%     UseCompass     (logical) — true if this vehicle uses compass assignment
%     CompassSector  (string)  — assigned compass sector ("NW","NE","SE","SW")
%                                empty string for non-compass vehicles
%
% COMPASS SECTOR MAPPING
% ----------------------
%   Effective corner  →  Compass sector  →  Carrier band [GHz]
%   FL                →  NE              →  78 – 79
%   FR                →  SE              →  79 – 80
%   RL                →  NW              →  77 – 78
%   RR                →  SW              →  80 – 81
%
%   For vehicles in opposing lanes (Invert = true in CornerInvertMap),
%   corner labels are rotated 180° before sector lookup:
%   FL↔RR,  FR↔RL.
%
% NOTES
% -----
%   - compassShare = 1 → all vehicles use compass (full COMPASS mode).
%   - compassShare = 0 → only the victim uses compass.
%   - Compass assignment pins Fc to the exact midpoint of the sector band.
%     F0 is then derived from Fc accounting for B, slope direction, and
%     the stepped frequency span.
%   - Non-compass vehicles use choose_f0_within_band (uniform draw),
%     identical to EXTENDED mode behaviour.
%
% SEE ALSO
%   ARIFOMA_generate_corner_parameters,
%   ARIFOMA_generate_fixed_pool_parameters,
%   ARIFOMA_build_parameter_table

% =========================================================================

%% ── RNG Initialisation ───────────────────────────────────────────────────
if cfg.UseFixedSeed
    rng(cfg.Seed);
end

veh     = vehicleIds(:);
corners = string(cfg.corners(:));
nV      = numel(veh);
nC      = numel(corners);
PRI_s   = cfg.CPI_duration_s;
nRows   = nV * nC;

%% ── Compass Vehicle Selection ────────────────────────────────────────────
%   Determine once per vehicle whether it uses compass frequency assignment.
%   The victim vehicle always uses compass; a random subset of the
%   remaining population is selected according to cfg.compassShare.

vehUsesCompass = select_compass_vehicles(veh, cfg);

%% ── Output Struct Pre-allocation ─────────────────────────────────────────
S = struct( ...
    'VehID',[], 'Corner',[], 'SensorID',[], 'HasSensor',[], ...
    'Scenario',[], 'Scheme',[], 'Waveform',[], 'WaveformClean',[], 'Dir',[], ...
    'UseCompass',[], 'CompassSector',[], ...
    'B',[], 'Tch',[], 'Nch',[], 'ChirpsPerStep',[], 'FreqIncrHz',[], ...
    'F0',[], 'Fc',[], ...
    'T_flyback_default',[], 'T_flyback_stepped',[], 'T_returnfly',[], ...
    'bandMin',[], 'bandMax',[], ...
    'Nseg',[], 'PermPattern',[]);

S   = repmat(S, nRows, 1);
row = 0;

%% ── Main Generation Loop ─────────────────────────────────────────────────
%   Outer loop — vehicles.
%   Inner loop — corners.
%   HasSensor and useCompassVeh are drawn once per vehicle.

for iV = 1:nV

    VehID         = veh(iV);
    vehHasSensors = (rand <= cfg.pHasSensor);
    useCompassVeh = vehUsesCompass(iV);

    for iC = 1:nC
        row    = row + 1;
        Corner = corners(iC);

        % ── Sensor ID Resolution ──────────────────────────────────────────
        SensorID = VehID * 100 + corner_index(Corner);

        if ~isempty(cfg.SensorIdMap)
            hit = cfg.SensorIdMap( ...
                cfg.SensorIdMap.VehID   == VehID & ...
                cfg.SensorIdMap.Corner  == Corner, :);
            if ~isempty(hit)
                SensorID = double(hit.SensorID(1));
            end
        end

        % ── Radar Class (Scenario) Selection ─────────────────────────────
        Scenario = pick_weighted( ...
            cfg.scenarioDist.(Corner).list, ...
            cfg.scenarioDist.(Corner).weights);

        sc = cfg.bandByScenario.(Scenario);

        bandMin = sc.bandMin;
        bandMax = sc.bandMax;
        Brange  = sc.Brange;
        Trange  = sc.Trange;

        % ── Sentinel Defaults (overwritten below if HasSensor = true) ─────
        Scheme = ""; Waveform = ""; WaveformClean = "";
        Dir = NaN;  B = NaN;  Tch = NaN;  Nch = NaN;
        ChirpsPerStep = NaN;  FreqIncrHz = NaN;
        F0 = NaN;   Fc = NaN;  CompassSector = "";

        Tfly_def  = min(sc.T_flyback_default);
        Tfly_step = cfg.FlybackSteppedMin_s;
        Tret      = cfg.ReturnFlyMin_s;

        if ~vehHasSensors && ~cfg.KeepIdForAbsent
            SensorID = NaN;
        end

        % ── Parameter Generation (sensor-equipped vehicles only) ──────────
        if vehHasSensors

            Scheme   = pick_weighted(cfg.SchemeList,   cfg.SchemeWeights);
            Waveform = pick_weighted(cfg.WaveformList, cfg.WaveformWeights);
            [WaveformClean, Dir] = normalize_Waveform(Waveform);

            % ── Sweep Bandwidth and Chirp Duration ────────────────────────
            B   = Brange(1) + rand * (Brange(2) - Brange(1));   % [Hz]
            Tch = Trange(1) + rand * (Trange(2) - Trange(1));   % [s]

            % ── Flyback Duration Sampling ─────────────────────────────────
            [lo, hi]  = order_pair(cfg.FlybackDefaultMin_s, cfg.FlybackDefaultMax_s);
            Tfly_def  = pick_real(lo, hi);

            [lo, hi]  = order_pair(cfg.FlybackSteppedMin_s, cfg.FlybackSteppedMax_s);
            Tfly_step = pick_real(lo, hi);

            [lo, hi]  = order_pair(cfg.ReturnFlyMin_s, cfg.ReturnFlyMax_s);
            Tret      = pick_real(lo, hi);

            % ── Number of Chirps per CPI ──────────────────────────────────
            if WaveformClean == "STEPPED FMCW"
                fly_s = Tfly_step;
            else
                fly_s = Tfly_def;
            end

            Nch = max(1, floor(PRI_s / max(eps, Tch + fly_s)));

            % ── Stepped FMCW Handling ─────────────────────────────────────
            if WaveformClean == "STEPPED FMCW"

                [lo, hi] = order_pair(cfg.StepCountMin, cfg.StepCountMax);
                lo       = max(1, floor(lo));
                hi       = max(lo, floor(hi));
                Nsteps   = pick_int(lo, hi);
                Nsteps   = min(Nsteps, Nch);

                [lo, hi] = order_pair(cfg.StepHzMin, cfg.StepHzMax);
                stepHz   = pick_real(lo, hi);

                [B, stepHz]   = fit_span_to_band(bandMin, bandMax, B, Nsteps, stepHz);
                FreqIncrHz    = stepHz;
                ChirpsPerStep = Nsteps;
                nStepForF0    = Nsteps;

            else
                ChirpsPerStep = NaN;
                FreqIncrHz    = 0;
                nStepForF0    = 1;
            end

            % ── Carrier Frequency Assignment ──────────────────────────────
            %   Compass vehicles: Fc is pinned to the exact midpoint of the
            %   sector band; F0 is derived from Fc via the sweep geometry.
            %   Non-compass vehicles: F0 drawn uniformly within the RF band.

            if useCompassVeh

                % Resolve effective compass sector for this (vehicle, corner).
                CompassSector = compass_sector_from_invert_map( ...
                    VehID, Corner, cfg.CornerInvertMap);

                % Target Fc = exact midpoint of the assigned sector band.
                fcBand    = cfg.compassFcBands.(CompassSector);
                Fc_target = 0.5 * (fcBand(1) + fcBand(2));

                [F0, Fc] = choose_f0_from_compass_center( ...
                    Fc_target, bandMin, bandMax, B, nStepForF0, Dir, FreqIncrHz);

            else
                % Random assignment — identical to EXTENDED mode behaviour.
                CompassSector = "";
                [F0, Fc] = choose_f0_within_band( ...
                    bandMin, bandMax, B, nStepForF0, Dir, FreqIncrHz);
            end

        end   % vehHasSensors

        % ── Segmentation Parameters (BlueFMCW) ────────────────────────────
        if isfinite(SensorID)
            Nseg        = 4 + mod(SensorID, 5);
            rs          = RandStream('mt19937ar', ...
                              'Seed', uint32(2654435761 * SensorID + iC));
            permPattern = randperm(rs, Nseg);
        else
            Nseg        = 4;
            permPattern = 1:Nseg;
        end

        % ── Store Row ─────────────────────────────────────────────────────
        S(row) = struct( ...
            'VehID',            VehID, ...
            'Corner',           Corner, ...
            'SensorID',         SensorID, ...
            'HasSensor',        vehHasSensors, ...
            'Scenario',         Scenario, ...
            'Scheme',           Scheme, ...
            'Waveform',         Waveform, ...
            'WaveformClean',    WaveformClean, ...
            'Dir',              Dir, ...
            'UseCompass',       useCompassVeh, ...
            'CompassSector',    CompassSector, ...
            'B',                B, ...
            'Tch',              Tch, ...
            'Nch',              Nch, ...
            'ChirpsPerStep',    ChirpsPerStep, ...
            'FreqIncrHz',       FreqIncrHz, ...
            'F0',               F0, ...
            'Fc',               Fc, ...
            'T_flyback_default', Tfly_def, ...
            'T_flyback_stepped', Tfly_step, ...
            'T_returnfly',      Tret, ...
            'bandMin',          bandMin, ...
            'bandMax',          bandMax, ...
            'Nseg',             Nseg, ...
            'PermPattern',      permPattern);

    end   % iC (corners)
end   % iV (vehicles)

Param = struct2table(S);

end   % ARIFOMA_generate_corner_parameters_compass

% =========================================================================
%% ── Local Helper Functions ───────────────────────────────────────────────
% =========================================================================

% -------------------------------------------------------------------------
function useCompassVeh = select_compass_vehicles(vehicleIds, cfg)
% SELECT_COMPASS_VEHICLES
%   Determine which vehicles receive compass-based frequency assignment.
%
%   The victim vehicle (cfg.VehID) is always assigned compass.
%   From the remaining population, exactly round(cfg.compassShare * N)
%   vehicles are randomly selected to also use compass.
%
%   INPUTS
%     vehicleIds  (double array) — all vehicle IDs in the scene
%     cfg         (struct)       — must contain cfg.VehID and cfg.compassShare
%
%   OUTPUT
%     useCompassVeh  (logical array, nV×1) — true for compass-assigned vehicles

    nV            = numel(vehicleIds);
    useCompassVeh = false(nV, 1);

    % The victim always uses compass.
    victimMask = false(nV, 1);
    if isfield(cfg, 'VehID') && ~isempty(cfg.VehID)
        victimMask            = (vehicleIds == cfg.VehID);
        useCompassVeh(victimMask) = true;
    end

    % Select a random subset of the non-victim population.
    otherIdx      = find(~victimMask);
    nOther        = numel(otherIdx);
    share         = max(0, min(1, cfg.compassShare));
    nCompassOther = round(share * nOther);

    if nCompassOther > 0
        pick              = otherIdx(randperm(nOther, nCompassOther));
        useCompassVeh(pick) = true;
    end
end

% -------------------------------------------------------------------------
function sector = compass_sector_from_invert_map(VehID, Corner, CornerInvertMap)
% COMPASS_SECTOR_FROM_INVERT_MAP
%   Resolve the compass sector for a given (vehicle, corner) pair.
%
%   Looks up the vehicle's inversion flag in CornerInvertMap and applies a
%   180-degree corner rotation if the vehicle is in the opposing lane.
%   The effective corner is then mapped to its compass sector.
%
%   Corner → Sector mapping (after inversion, if applicable):
%     FL → NE  |  FR → SE  |  RL → NW  |  RR → SW
%
%   INPUTS
%     VehID          (double) — vehicle ID to look up
%     Corner         (string) — original corner label
%     CornerInvertMap (table) — {VehID, Invert} produced by
%                               ARIFOMA_build_parameter_table
%
%   OUTPUT
%     sector  (string) — compass sector: "NW", "NE", "SE", or "SW"

    invFlag    = lookup_corner_invert_flag(VehID, CornerInvertMap);
    CornerEff  = apply_corner_inversion_180(Corner, invFlag);

    switch string(CornerEff)
        case "FL",    sector = "NE";
        case "FR",    sector = "SE";
        case "RL",    sector = "NW";
        case "RR",    sector = "SW";
        otherwise
            error('ARIFOMA:compassSector:unknownCorner', ...
                  'Unknown corner label: "%s".', string(Corner));
    end
end

% -------------------------------------------------------------------------
function invFlag = lookup_corner_invert_flag(VehID, CornerInvertMap)
% LOOKUP_CORNER_INVERT_FLAG
%   Return the inversion flag for a vehicle from CornerInvertMap.
%
%   Returns false if the map is empty or the vehicle is not found.
%
%   INPUTS
%     VehID           (double) — vehicle ID to look up
%     CornerInvertMap (table)  — {VehID, Invert}
%
%   OUTPUT
%     invFlag  (logical) — true if this vehicle's corners should be inverted

    invFlag = false;

    if isempty(CornerInvertMap)
        return;
    end

    m = (CornerInvertMap.VehID == VehID);
    if any(m)
        invFlag = logical(CornerInvertMap.Invert(find(m, 1, 'first')));
    end
end

% -------------------------------------------------------------------------
function CornerOut = apply_corner_inversion_180(CornerIn, invFlag)
% APPLY_CORNER_INVERSION_180
%   Rotate corner labels 180 degrees for vehicles in opposing traffic lanes.
%
%   When a vehicle travels in the opposing direction, its physical sensor
%   corners are geometrically mirrored relative to the ego vehicle's frame.
%   The rotation map is:  FL ↔ RR,  FR ↔ RL.
%
%   INPUTS
%     CornerIn  (string)  — original corner label
%     invFlag   (logical) — true if inversion should be applied
%
%   OUTPUT
%     CornerOut  (string) — effective corner label after inversion

    CornerIn = string(CornerIn);

    if ~invFlag
        CornerOut = CornerIn;
        return;
    end

    switch CornerIn
        case "FL",    CornerOut = "RR";
        case "FR",    CornerOut = "RL";
        case "RL",    CornerOut = "FR";
        case "RR",    CornerOut = "FL";
        otherwise
            error('ARIFOMA:cornerInversion:unknownCorner', ...
                  'Unknown corner label: "%s".', CornerIn);
    end
end

% -------------------------------------------------------------------------
function [F0, Fc] = choose_f0_from_compass_center( ...
        Fc_target, bandMin, bandMax, B, Nsteps, dir, stepHz)
% CHOOSE_F0_FROM_COMPASS_CENTER
%   Derive F0 by centring the sweep span on a target compass carrier frequency.
%
%   The sweep is placed so that its centre frequency matches Fc_target as
%   closely as possible.  If the resulting base frequency falls outside the
%   RF band boundaries, it is clipped to the nearest valid position.
%
%   INPUTS
%     Fc_target          (double) — target centre frequency [Hz]
%     bandMin, bandMax   (double) — RF band limits [Hz]
%     B                  (double) — sweep bandwidth [Hz]
%     Nsteps             (double) — number of stepped sub-sweeps (1 for standard)
%     dir                (double) — chirp slope sign (+1 or -1)
%     stepHz             (double) — frequency step size [Hz] (0 for standard FMCW)
%
%   OUTPUTS
%     F0  (double) — starting carrier frequency [Hz]
%     Fc  (double) — achieved centre frequency [Hz]

    if ~isfinite(stepHz) || stepHz < 0,  stepHz = 0;  end
    if ~isfinite(Nsteps) || Nsteps < 1,  Nsteps = 1;  end

    span    = max(B + (Nsteps - 1) * stepHz, B);
    bwAvail = bandMax - bandMin;

    % Safety clip on span.
    if span > bwAvail
        span = bwAvail;
    end

    % Centre the span on Fc_target, then clip base to valid range.
    base = Fc_target - span / 2;
    base = min(max(base, bandMin), bandMax - span);

    if dir >= 0
        F0 = base;
    else
        F0 = base + B;
    end

    Fc = base + span / 2;
end

% -------------------------------------------------------------------------
function [B_adj, stepHz_adj] = fit_span_to_band(bandMin, bandMax, B, Nch, stepHz)
% FIT_SPAN_TO_BAND  Clip (B, stepHz) so the total stepped span fits the RF band.
%
%   Enforces:  span = B + (Nch-1)*stepHz  <=  (bandMax - bandMin)
%
%   Strategy:
%     1. Clip B to bwAvail if it alone exceeds the available bandwidth.
%     2. Derive the maximum allowable stepHz given the clipped B and Nch.
%     3. Clip stepHz to this maximum.
%     4. Apply a final B clip as a numerical safety net.
%
%   INPUTS
%     bandMin, bandMax  (double) — RF band limits [Hz]
%     B                (double) — sweep bandwidth [Hz]
%     Nch              (double) — number of stepped sub-sweeps
%     stepHz           (double) — frequency step size [Hz]
%
%   OUTPUTS
%     B_adj      (double) — clipped sweep bandwidth [Hz]
%     stepHz_adj (double) — clipped frequency step size [Hz]

    bandMin = double(bandMin);  bandMax = double(bandMax);
    B       = double(B);        Nch     = double(Nch);
    stepHz  = double(stepHz);

    bwAvail = bandMax - bandMin;

    if Nch <= 1
        B_adj      = min(B, bwAvail);
        stepHz_adj = 0;
        return;
    end

    if B > bwAvail,  B = bwAvail;  end

    maxStepHz  = max(0, (bwAvail - B) / (Nch - 1));
    stepHz_adj = min(max(0, stepHz), maxStepHz);

    span = B + (Nch - 1) * stepHz_adj;
    if span > bwAvail
        B = max(0, bwAvail - (Nch - 1) * stepHz_adj);
    end

    B_adj = B;
end

% -------------------------------------------------------------------------
function [wfClean, dir] = normalize_Waveform(label, dirExisting)
% NORMALIZE_WAVEFORM  Map a waveform label string to canonical type and slope sign.
%
%   INPUTS
%     label       (string) — raw waveform label from cfg.WaveformList
%     dirExisting (double) — optional existing slope sign to preserve
%
%   OUTPUTS
%     wfClean  (string) — "FMCW" or "STEPPED FMCW"
%     dir      (double) — chirp slope sign: +1 or -1

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

% -------------------------------------------------------------------------
function x = pick_weighted(values, weights)
% PICK_WEIGHTED  Draw one element from `values` according to `weights`.

    values = string(values(:));
    w      = weights(:) / sum(weights(:));
    x      = values(find(rand <= cumsum(w), 1, 'first'));
end

% -------------------------------------------------------------------------
function [F0, Fc] = choose_f0_within_band(bandMin, bandMax, B, Nch, dir, stepHz)
% CHOOSE_F0_WITHIN_BAND  Draw a starting carrier frequency uniformly within the RF band.
%
%   Used for non-compass vehicles — behaviour is identical to EXTENDED mode.
%
%   INPUTS
%     bandMin, bandMax  (double) — RF band limits [Hz]
%     B                (double) — sweep bandwidth [Hz]
%     Nch              (double) — number of stepped sub-sweeps
%     dir              (double) — chirp slope sign (+1 or -1)
%     stepHz           (double) — frequency step size [Hz]
%
%   OUTPUTS
%     F0  (double) — starting carrier frequency [Hz]
%     Fc  (double) — centre carrier frequency [Hz]

    if ~isfinite(stepHz) || stepHz <= 0,  stepHz = 0;  end

    span = max(B + (Nch - 1) * stepHz, B);

    if span > (bandMax - bandMin)
        B    = min(B, (bandMax - bandMin) - max(0, (Nch - 1) * stepHz));
        span = B + (Nch - 1) * stepHz;
    end

    base = bandMin + rand * max(0, (bandMax - bandMin) - span);

    if dir >= 0
        F0 = base;
    else
        F0 = base + B;
    end

    Fc = base + B / 2;
end

% -------------------------------------------------------------------------
function [lo, hi] = order_pair(a, b)
% ORDER_PAIR  Return (lo, hi) regardless of input order.

    lo = min(double(a), double(b));
    hi = max(double(a), double(b));
end

% -------------------------------------------------------------------------
function x = pick_real(lo, hi)
% PICK_REAL  Draw a scalar uniformly from [lo, hi].

    if lo == hi,  x = lo;
    else,         x = lo + rand * (hi - lo);
    end
end

% -------------------------------------------------------------------------
function k = pick_int(lo, hi)
% PICK_INT  Draw an integer uniformly from {lo, lo+1, ..., hi}.

    if lo == hi,  k = lo;
    else,         k = randi([lo hi]);
    end
end

% -------------------------------------------------------------------------
function idx = corner_index(Corner)
% CORNER_INDEX  Map a corner label string to its integer index.
%
%   FL → 1,  FR → 2,  RL → 3,  RR → 4,  F → 5
%   Returns 0 for unrecognised labels.

    map = struct('FL',1, 'FR',2, 'RL',3, 'RR',4, 'F',5);
    if isfield(map, char(Corner))
        idx = map.(char(Corner));
    else
        idx = 0;
    end
end