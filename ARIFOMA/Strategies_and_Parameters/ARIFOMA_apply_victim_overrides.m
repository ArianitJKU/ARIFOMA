function Parameters = ARIFOMA_apply_victim_overrides(Parameters, VehID, cfg)
% ARIFOMA_APPLY_VICTIM_OVERRIDES
%   Apply victim-specific scenario, waveform, and scheme overrides and
%   regenerate all dependent radar parameters for the victim vehicle only.
%
% PURPOSE
% -------
%   During simulation setup, the victim vehicle may need to be assigned
%   a specific radar class (scenario), waveform type, or mitigation scheme
%   that differs from the randomly drawn population parameters.  This
%   function enforces those overrides and re-derives all downstream
%   parameters (bandwidth, chirp duration, flyback, Nch, stepped FMCW
%   parameters, F0, Fc) so that the parameter table remains internally
%   consistent after the override.
%
%   The function modifies only the rows of Parameters corresponding to the
%   victim vehicle (VehID) that have HasSensor = true.  All other rows are
%   left unchanged.
%
% SIX-STEP OVERRIDE PIPELINE (per victim corner row)
% ---------------------------------------------------
%   1. Scenario override    — apply cfg.victimscenarioPerCorner per corner
%                             and update bandMin / bandMax accordingly
%   2. Waveform / scheme    — overwrite from cfg.adaptWaveformVictim and
%                             cfg.adaptSchemeVictim; renormalise WaveformClean
%                             and Dir
%   3. Analog parameters    — redraw B and Tch from the new scenario ranges;
%                             redraw flyback / return-fly timing
%   4. Nch recomputation    — derive number of chirps per CPI from the
%                             updated Tch and flyback duration
%   5. Stepped FMCW params  — redraw Nsteps and stepHz; clip to RF band;
%                             update ChirpsPerStep and FreqIncrHz
%   6. F0 / Fc              — draw carrier frequency within the new scenario
%                             band using the updated sweep geometry
%
% INPUTS
%   Parameters  (table)  — full per-sensor parameter table produced by
%                          ARIFOMA_build_parameter_table; modified in-place
%   VehID       (double) — victim vehicle ID; only rows with
%                          Parameters.VehID == VehID and HasSensor == true
%                          are modified
%   cfg         (struct) — global configuration struct; relevant fields:
%
%     Scenario override
%       cfg.victimscenarioPerCorner  (struct) — per-corner scenario strings,
%                                               fields FL, FR, RL, RR;
%                                               if absent, current scenario
%                                               is preserved
%       cfg.bandByScenario           (struct) — per-scenario RF and waveform
%                                               parameter ranges
%
%     Waveform / scheme override
%       cfg.adaptSchemeVictim        (string) — scheme to force onto victim;
%                                               applied if non-empty
%       cfg.adaptWaveformVictim      (string) — waveform type to force onto
%                                               victim; applied if non-empty
%
%     Timing
%       cfg.CPI_duration_s           (double) — CPI / PRI duration [s]
%       cfg.FlybackDefaultMin/Max_s  (double) — default flyback bounds [s]
%       cfg.FlybackSteppedMin/Max_s  (double) — stepped flyback bounds [s]
%       cfg.ReturnFlyMin/Max_s       (double) — return-fly bounds [s]
%
%     Stepped FMCW
%       cfg.StepCountMin/Max         (double) — step count bounds
%       cfg.StepHzMin/Max            (double) — step size bounds [Hz]
%
% OUTPUT
%   Parameters  (table) — parameter table with victim rows updated;
%                         all six parameter groups are internally consistent
%
% NOTES
% -----
%   - If VehID is not found in Parameters (no sensor-equipped rows), the
%     function returns immediately without modification.
%   - This function is only called when cfg.paramMode ~= "COMPASS".
%     In COMPASS mode, carrier frequencies are governed by directional
%     band assignment and overrides are applied differently.
%   - Analog parameters (B, Tch, flyback) are redrawn randomly from the
%     new scenario ranges.  To fix these values, set the range bounds
%     equal (e.g. Brange = [B_fixed, B_fixed]).
%
% SEE ALSO
%   ARIFOMA_build_parameter_table,
%   ARIFOMA_generate_corner_parameters,
%   main_sim

% =========================================================================

%% ── Identify Victim Rows ─────────────────────────────────────────────────
%   Only rows belonging to the victim vehicle with an active sensor are
%   processed.  All other rows are left untouched.

isVict  = (Parameters.VehID == VehID) & (Parameters.HasSensor);
idxVict = find(isVict);

if isempty(idxVict)
    return;   % victim not found or has no sensors — nothing to override
end

%% ── Per-Corner Override Loop ─────────────────────────────────────────────
for k = idxVict(:).'

    corner = string(Parameters.Corner(k));

    % =====================================================================
    % Step 1: Scenario Override Per Corner
    % =====================================================================
    %   Apply a corner-specific radar class (MRR / SRR / LRR) from
    %   cfg.victimscenarioPerCorner if the field exists for this corner.
    %   If the field is absent, the current scenario assignment is preserved.

    if isfield(cfg, 'victimscenarioPerCorner') && ...
            isfield(cfg.victimscenarioPerCorner, corner)
        scenarioVict = string(cfg.victimscenarioPerCorner.(corner));
    else
        scenarioVict = string(Parameters.Scenario(k));
    end

    Parameters.Scenario(k) = scenarioVict;

    % Update RF band limits from the new scenario definition.
    sc = cfg.bandByScenario.(scenarioVict);
    Parameters.bandMin(k) = sc.bandMin;
    Parameters.bandMax(k) = sc.bandMax;

    % =====================================================================
    % Step 2: Waveform and Scheme Override
    % =====================================================================
    %   Overwrite scheme and waveform from cfg if the override fields are
    %   non-empty.  WaveformClean and Dir are re-derived from the new
    %   waveform label to keep the table consistent.

    if isfield(cfg, 'adaptSchemeVictim') && ~isempty(cfg.adaptSchemeVictim)
        Parameters.Scheme(k) = string(cfg.adaptSchemeVictim);
    end

    if isfield(cfg, 'adaptWaveformVictim') && ~isempty(cfg.adaptWaveformVictim)
        Parameters.Waveform(k) = string(cfg.adaptWaveformVictim);
    end

    % Re-derive canonical waveform type and slope direction from the label.
    [wfClean, dirVal]        = normalize_Waveform(Parameters.Waveform(k));
    Parameters.WaveformClean(k) = wfClean;
    Parameters.Dir(k)           = dirVal;

    % =====================================================================
    % Step 3: Regenerate Scenario-Dependent Analog Parameters
    % =====================================================================
    %   Redraw bandwidth and chirp duration from the new scenario ranges.
    %   Flyback and return-fly durations are also redrawn to maintain
    %   consistency with the generator functions.

    Brange = sc.Brange;
    Trange = sc.Trange;

    Parameters.B(k)   = Brange(1) + rand * (Brange(2) - Brange(1));   % [Hz]
    Parameters.Tch(k) = Trange(1) + rand * (Trange(2) - Trange(1));   % [s]

    [lo, hi] = order_pair(cfg.FlybackDefaultMin_s, cfg.FlybackDefaultMax_s);
    Parameters.T_flyback_default(k) = pick_real(lo, hi);

    [lo, hi] = order_pair(cfg.FlybackSteppedMin_s, cfg.FlybackSteppedMax_s);
    Parameters.T_flyback_stepped(k) = pick_real(lo, hi);

    [lo, hi] = order_pair(cfg.ReturnFlyMin_s, cfg.ReturnFlyMax_s);
    Parameters.T_returnfly(k) = pick_real(lo, hi);

    % =====================================================================
    % Step 4: Recompute Number of Chirps per CPI (Nch)
    % =====================================================================
    %   Nch is derived from the CPI duration and the effective PRI length
    %   (Tch + flyback).  The flyback duration used depends on whether the
    %   waveform is stepped or standard.

    if Parameters.WaveformClean(k) == "STEPPED FMCW"
        fly_s = Parameters.T_flyback_stepped(k);
    else
        fly_s = Parameters.T_flyback_default(k);
    end

    Parameters.Nch(k) = max(1, floor( ...
        cfg.CPI_duration_s / max(eps, Parameters.Tch(k) + fly_s)));

    % =====================================================================
    % Step 5: Recompute Stepped FMCW Parameters
    % =====================================================================
    %   For stepped waveforms, redraw the step count and step size and clip
    %   the resulting swept span to fit within the RF band.
    %   For standard FMCW, clear the stepped fields.

    if Parameters.WaveformClean(k) == "STEPPED FMCW"

        % Robust step count sampling with floor guards.
        [lo, hi] = order_pair(cfg.StepCountMin, cfg.StepCountMax);
        lo       = max(1, floor(lo));
        hi       = max(lo, floor(hi));
        Nsteps   = pick_int(lo, hi);
        Nsteps   = min(Nsteps, Parameters.Nch(k));   % cannot exceed chirp count

        % Robust step size sampling.
        [lo, hi] = order_pair(cfg.StepHzMin, cfg.StepHzMax);
        stepHz   = pick_real(lo, hi);

        % Clip (B, stepHz) so the total swept span fits the RF band.
        [Bnew, stepHz] = fit_span_to_band( ...
            Parameters.bandMin(k), Parameters.bandMax(k), ...
            Parameters.B(k), Nsteps, stepHz);

        Parameters.B(k)             = Bnew;
        Parameters.ChirpsPerStep(k) = Nsteps;
        Parameters.FreqIncrHz(k)    = stepHz;
        nStepForF0                  = Nsteps;

    else
        % Standard FMCW — clear stepped parameters.
        Parameters.ChirpsPerStep(k) = NaN;
        Parameters.FreqIncrHz(k)    = 0;
        nStepForF0                  = 1;
    end

    % =====================================================================
    % Step 6: Recompute F0 and Fc Within the New Scenario Band
    % =====================================================================
    %   Draw a carrier frequency uniformly within the available margin of
    %   the new scenario RF band, accounting for the full swept span.

    [Parameters.F0(k), Parameters.Fc(k)] = choose_f0_within_band( ...
        Parameters.bandMin(k), ...
        Parameters.bandMax(k), ...
        Parameters.B(k), ...
        nStepForF0, ...
        Parameters.Dir(k), ...
        Parameters.FreqIncrHz(k));

end   % per-corner override loop

end   % ARIFOMA_apply_victim_overrides

% =========================================================================
%% ── Local Helper Functions ───────────────────────────────────────────────
% =========================================================================

% -------------------------------------------------------------------------
function [wfClean, dir] = normalize_Waveform(label, dirExisting)
% NORMALIZE_WAVEFORM  Map a waveform label string to canonical type and slope sign.
%
%   INPUTS
%     label       (string) — raw waveform label from Parameters.Waveform
%     dirExisting (double) — optional existing slope sign to preserve (+1 or -1)
%
%   OUTPUTS
%     wfClean  (string) — canonical type: "FMCW" or "STEPPED FMCW"
%     dir      (double) — chirp slope sign: +1 (positive) or -1 (negative)

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
function [lo, hi] = order_pair(a, b)
% ORDER_PAIR  Return (lo, hi) regardless of input order.
%
%   Guards against swapped min/max entries in cfg.

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
function [F0, Fc] = choose_f0_within_band(bandMin, bandMax, B, Nch, dir, stepHz)
% CHOOSE_F0_WITHIN_BAND  Draw a starting carrier frequency within the RF band.
%
%   F0 is drawn uniformly from the range [bandMin, bandMax - span], where
%   span = B + (Nch-1)*stepHz accounts for the full stepped frequency
%   excursion.  The sign of dir determines whether F0 anchors the low
%   (positive slope) or high (negative slope) end of the sweep.
%
%   INPUTS
%     bandMin, bandMax  (double) — RF band limits [Hz]
%     B                (double) — sweep bandwidth [Hz]
%     Nch              (double) — number of stepped sub-sweeps (1 for standard)
%     dir              (double) — chirp slope sign (+1 or -1)
%     stepHz           (double) — frequency step size [Hz]
%
%   OUTPUTS
%     F0  (double) — starting carrier frequency [Hz]
%     Fc  (double) — centre carrier frequency [Hz]

    if ~isfinite(stepHz) || stepHz <= 0,  stepHz = 0;  end

    span = max(B + (Nch - 1) * stepHz, B);

    % Safety clip — should already be handled by fit_span_to_band.
    if span > (bandMax - bandMin)
        B    = min(B, (bandMax - bandMin) - max(0, (Nch - 1) * stepHz));
        span = B + (Nch - 1) * stepHz;
    end

    margin = max(0, (bandMax - bandMin) - span);
    base   = bandMin + rand * margin;

    if dir >= 0
        F0 = base;       % Positive slope: sweep starts at base
    else
        F0 = base + B;   % Negative slope: sweep starts at base + B
    end

    Fc = base + B / 2;   % Centre frequency (slope-sign independent)
end