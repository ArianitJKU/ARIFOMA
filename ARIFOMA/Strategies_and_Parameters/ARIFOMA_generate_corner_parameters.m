function Param = ARIFOMA_generate_corner_parameters(vehicleIds, cfg)
% ARIFOMA_GENERATE_CORNER_PARAMETERS
%   Per-corner FMCW radar parameter generator (no compass logic).
%
% PURPOSE
% -------
%   Assigns independent, randomly drawn FMCW waveform parameters to every
%   (vehicle, corner) sensor instance in the scene.  Unlike the COMPASS
%   mode generator, carrier frequencies are drawn freely from the full RF
%   band without directional sub-band constraints.
%
% FEATURES
% --------
%   - Number of chirps per CPI (Nch) derived directly from PRI duration
%   - Robust to swapped min/max configuration inputs (order_pair guards)
%   - Stepped FMCW support with band-clipped span enforcement
%   - Compatible with GUI-driven parameter configuration
%   - Sensor ID mapping via optional cfg.SensorIdMap lookup table
%
% INPUTS
%   vehicleIds  (double array) — column vector of unique vehicle IDs
%   cfg         (struct)       — global configuration struct; required fields:
%
%     RNG
%       cfg.UseFixedSeed       (logical) — seed the RNG before generation
%       cfg.Seed               (double)  — RNG seed value
%
%     Geometry
%       cfg.corners            (string array) — corner labels, e.g. ["FL","FR","RL","RR"]
%       cfg.SensorIdMap        (table)   — optional external sensor-ID lookup;
%                                          columns: VehID, Corner, SensorID
%       cfg.KeepIdForAbsent    (logical) — retain sensor IDs for sensor-less corners
%       cfg.pHasSensor         (double)  — probability a vehicle carries sensors [0,1]
%
%     Timing
%       cfg.CPI_duration_s     (double)  — full PRI / CPI duration [s]
%       cfg.FlybackDefaultMin_s (double) — default flyback lower bound [s]
%       cfg.FlybackDefaultMax_s (double) — default flyback upper bound [s]
%       cfg.FlybackSteppedMin_s (double) — stepped flyback lower bound [s]
%       cfg.FlybackSteppedMax_s (double) — stepped flyback upper bound [s]
%       cfg.ReturnFlyMin_s     (double)  — return-fly lower bound [s]
%       cfg.ReturnFlyMax_s     (double)  — return-fly upper bound [s]
%
%     RF Band and Waveform
%       cfg.scenarioDist       (struct)  — per-corner radar-class distributions
%       cfg.bandByScenario     (struct)  — per-class RF and waveform parameter ranges
%       cfg.SchemeList         (string array) — mitigation scheme pool
%       cfg.SchemeWeights      (double array) — scheme selection weights
%       cfg.WaveformList       (string array) — waveform type pool
%       cfg.WaveformWeights    (double array) — waveform selection weights
%
%     Stepped FMCW
%       cfg.StepCountMin       (double)  — minimum number of sub-sweep steps
%       cfg.StepCountMax       (double)  — maximum number of sub-sweep steps
%       cfg.StepHzMin          (double)  — minimum frequency step size [Hz]
%       cfg.StepHzMax          (double)  — maximum frequency step size [Hz]
%
% OUTPUT
%   Param  (table) — parameter table with one row per (vehicle, corner) pair.
%
%   Schema (one column per field):
%     VehID, Corner, SensorID, HasSensor
%     Scenario, Scheme, Waveform, WaveformClean, Dir
%     B, Tch, Nch, ChirpsPerStep, FreqIncrHz
%     F0, Fc
%     T_flyback_default, T_flyback_stepped, T_returnfly
%     bandMin, bandMax
%     Nseg, PermPattern
%
% NOTES
% -----
%   - Each vehicle independently draws a HasSensor flag from cfg.pHasSensor.
%     All corners of a sensor-less vehicle receive NaN waveform fields.
%   - Segmentation (Nseg, PermPattern) is computed deterministically from
%     SensorID so results are reproducible without a global RNG seed.
%   - For compass-based carrier-frequency assignment, use
%     ARIFOMA_generate_compass_parameters instead.
%
% SEE ALSO
%   ARIFOMA_generate_compass_parameters,
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

PRI_s = cfg.CPI_duration_s;   % Full PRI duration used to derive Nch [s]

%% ── Output Struct Pre-allocation ─────────────────────────────────────────
%   Pre-allocating a struct array and converting to table at the end is
%   significantly faster than growing a table row-by-row inside the loop.

nRows = nV * nC;

S = struct( ...
    'VehID',[], 'Corner',[], 'SensorID',[], 'HasSensor',[], ...
    'Scenario',[], 'Scheme',[], 'Waveform',[], 'WaveformClean',[], 'Dir',[], ...
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
%   HasSensor is drawn once per vehicle so all corners share the same flag.

for iV = 1:nV

    VehID         = veh(iV);
    vehHasSensors = (rand <= cfg.pHasSensor);   % one draw per vehicle

    for iC = 1:nC
        row    = row + 1;
        Corner = corners(iC);

        % ── Sensor ID Resolution ──────────────────────────────────────────
        %   Default: deterministic encoding as VehID*100 + corner index.
        %   Override: look up in cfg.SensorIdMap if provided.
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
        %   Drawn from the per-corner weighted distribution.
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
        F0 = NaN;   Fc = NaN;

        Tfly_def  = min(sc.T_flyback_default);
        Tfly_step = cfg.FlybackSteppedMin_s;
        Tret      = cfg.ReturnFlyMin_s;

        % Optionally clear sensor ID for sensor-less corners.
        if ~vehHasSensors && ~cfg.KeepIdForAbsent
            SensorID = NaN;
        end

        % ── Parameter Generation (sensor-equipped vehicles only) ──────────
        if vehHasSensors

            % Scheme and waveform type drawn from configured pools.
            Scheme   = pick_weighted(cfg.SchemeList,   cfg.SchemeWeights);
            Waveform = pick_weighted(cfg.WaveformList, cfg.WaveformWeights);
            [WaveformClean, Dir] = normalize_Waveform(Waveform);

            % ── Sweep Bandwidth and Chirp Duration ────────────────────────
            B   = Brange(1) + rand * (Brange(2) - Brange(1));   % [Hz]
            Tch = Trange(1) + rand * (Trange(2) - Trange(1));   % [s]

            % ── Flyback Duration Sampling ─────────────────────────────────
            %   order_pair() guards against swapped min/max in cfg.
            [lo, hi]  = order_pair(cfg.FlybackDefaultMin_s, cfg.FlybackDefaultMax_s);
            Tfly_def  = pick_real(lo, hi);

            [lo, hi]  = order_pair(cfg.FlybackSteppedMin_s, cfg.FlybackSteppedMax_s);
            Tfly_step = pick_real(lo, hi);

            [lo, hi]  = order_pair(cfg.ReturnFlyMin_s, cfg.ReturnFlyMax_s);
            Tret      = pick_real(lo, hi);

            % Select which flyback duration governs the PRI length.
            if WaveformClean == "STEPPED FMCW"
                fly_s = Tfly_step;
            else
                fly_s = Tfly_def;
            end

            % ── Number of Chirps per CPI ──────────────────────────────────
            %   Derived from PRI_s so the CPI is fully covered.
            Nch = max(1, floor(PRI_s / max(eps, Tch + fly_s)));

            % ── Stepped FMCW Handling ─────────────────────────────────────
            if WaveformClean == "STEPPED FMCW"

                % Robust step count sampling (floor guards non-integer cfg).
                [lo, hi] = order_pair(cfg.StepCountMin, cfg.StepCountMax);
                lo       = max(1, floor(lo));
                hi       = max(lo, floor(hi));
                Nsteps   = pick_int(lo, hi);
                Nsteps   = min(Nsteps, Nch);   % Cannot exceed chirp count

                % Robust step size sampling.
                [lo, hi] = order_pair(cfg.StepHzMin, cfg.StepHzMax);
                stepHz   = pick_real(lo, hi);

                % Clip (B, stepHz) so the total stepped span fits the RF band.
                [B, stepHz] = fit_span_to_band(bandMin, bandMax, B, Nsteps, stepHz);

                ChirpsPerStep = Nsteps;
                FreqIncrHz    = stepHz;

            else
                % Standard FMCW — no sub-sweep stepping.
                ChirpsPerStep = NaN;
                FreqIncrHz    = 0;
            end

            % ── Carrier Frequency Selection ───────────────────────────────
            %   F0 is drawn uniformly within the available margin inside the
            %   RF band, accounting for the full swept span.
            [F0, Fc] = choose_f0_within_band( ...
                bandMin, bandMax, B, ChirpsPerStep, Dir, FreqIncrHz);

        end   % vehHasSensors

        % ── Segmentation Parameters (BlueFMCW) ────────────────────────────
        %   Nseg and PermPattern are derived deterministically from SensorID
        %   so that BlueFMCW permutations are reproducible across runs
        %   without relying on the global RNG state.
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

end   % ARIFOMA_generate_corner_parameters

% =========================================================================
%% ── Local Helper Functions ───────────────────────────────────────────────
% =========================================================================

% -------------------------------------------------------------------------
function [B_adj, stepHz_adj] = fit_span_to_band(bandMin, bandMax, B, Nch, stepHz)
% FIT_SPAN_TO_BAND  Clip (B, stepHz) so the total stepped span fits the RF band.
%
%   Enforces:  span = B + (Nch-1)*stepHz  <=  (bandMax - bandMin)
%
%   Strategy:
%     1. Clip B to bwAvail if it already exceeds the available bandwidth.
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

    % Degenerate case: single sub-sweep, no span accumulation.
    if Nch <= 1
        B_adj      = min(B, bwAvail);
        stepHz_adj = 0;
        return;
    end

    % Step 1 — clip B if it alone exceeds the available bandwidth.
    if B > bwAvail
        B = bwAvail;
    end

    % Step 2 — derive and enforce the maximum allowable step size.
    maxStepHz  = max(0, (bwAvail - B) / (Nch - 1));
    stepHz_adj = min(max(0, stepHz), maxStepHz);

    % Step 3 — final safety clip on B for numerical edge cases.
    span = B + (Nch - 1) * stepHz_adj;
    if span > bwAvail
        B = max(0, bwAvail - (Nch - 1) * stepHz_adj);
    end

    B_adj = B;
end

% -------------------------------------------------------------------------
function [wfClean, dir] = normalize_Waveform(label, dirExisting)
% NORMALIZE_WAVEFORM  Map a waveform label string to a canonical type and slope sign.
%
%   INPUTS
%     label       (string) — raw waveform label from cfg.WaveformList
%     dirExisting (double) — optional existing slope sign to preserve (+1 or -1)
%
%   OUTPUTS
%     wfClean  (string) — canonical type: "FMCW" or "STEPPED FMCW"
%     dir      (double) — chirp slope sign: +1 (positive) or -1 (negative)

    L = upper(string(label));
    switch L
        case "POSITIVE-SLOPE"
            wfClean = "FMCW";         dir = +1;
        case "NEGATIVE-SLOPE"
            wfClean = "FMCW";         dir = -1;
        case "POSITIVE-SLOPE STEPPED"
            wfClean = "STEPPED FMCW"; dir = +1;
        case "NEGATIVE-SLOPE STEPPED"
            wfClean = "STEPPED FMCW"; dir = -1;
        otherwise
            wfClean = "FMCW";         dir = +1;
    end

    % Preserve an existing slope direction if explicitly provided.
    if nargin >= 2 && isfinite(dirExisting) && dirExisting ~= 0
        dir = sign(dirExisting);
    end
end

% -------------------------------------------------------------------------
function x = pick_weighted(values, weights)
% PICK_WEIGHTED  Draw one element from `values` according to `weights`.
%
%   Uses inverse CDF sampling.  Weights are normalised internally and do
%   not need to sum to 1.
%
%   INPUTS
%     values   (string array) — candidate values
%     weights  (double array) — non-negative selection weights
%
%   OUTPUT
%     x  (string) — selected element

    values = string(values(:));
    w      = weights(:) / sum(weights(:));
    cdf    = cumsum(w);
    x      = values(find(rand <= cdf, 1, 'first'));
end

% -------------------------------------------------------------------------
function [F0, Fc] = choose_f0_within_band(bandMin, bandMax, B, Nch, dir, stepHz)
% CHOOSE_F0_WITHIN_BAND  Draw a starting carrier frequency within the RF band.
%
%   F0 is drawn uniformly from the range [bandMin, bandMax - span], where
%   span = B + (Nch-1)*stepHz accounts for the full stepped frequency
%   excursion.  The sign of `dir` determines whether F0 anchors the low
%   (positive slope) or high (negative slope) end of the sweep.
%
%   INPUTS
%     bandMin, bandMax  (double) — RF band limits [Hz]
%     B                (double) — sweep bandwidth [Hz]
%     Nch              (double) — number of stepped sub-sweeps (NaN → 1)
%     dir              (double) — chirp slope sign (+1 or -1)
%     stepHz           (double) — frequency step size [Hz]
%
%   OUTPUTS
%     F0  (double) — starting carrier frequency [Hz]
%     Fc  (double) — centre carrier frequency [Hz]

    if ~isfinite(stepHz) || stepHz <= 0
        stepHz = 0;
    end

    span = max(B + (Nch - 1) * stepHz, B);

    % Safety clip — should already be handled by fit_span_to_band.
    if span > (bandMax - bandMin)
        B    = min(B, (bandMax - bandMin) - max(0, (Nch - 1) * stepHz));
        span = B + (Nch - 1) * stepHz;
    end

    % Uniform draw within the available margin.
    margin = max(0, (bandMax - bandMin) - span);
    base   = bandMin + rand * margin;

    if dir >= 0
        F0 = base;           % Positive slope: sweep starts at base
    else
        F0 = base + B;       % Negative slope: sweep starts at base + B
    end

    Fc = base + B / 2;       % Centre frequency (independent of slope sign)
end

% -------------------------------------------------------------------------
function [lo, hi] = order_pair(a, b)
% ORDER_PAIR  Return (lo, hi) regardless of input order.
%
%   Prevents silent errors when cfg min/max fields are accidentally swapped.

    lo = min(double(a), double(b));
    hi = max(double(a), double(b));
end

% -------------------------------------------------------------------------
function x = pick_real(lo, hi)
% PICK_REAL  Draw a scalar uniformly from [lo, hi].

    if lo == hi
        x = lo;
    else
        x = lo + rand * (hi - lo);
    end
end

% -------------------------------------------------------------------------
function k = pick_int(lo, hi)
% PICK_INT  Draw an integer uniformly from {lo, lo+1, ..., hi}.

    if lo == hi
        k = lo;
    else
        k = randi([lo hi]);
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