function [corner, radarCtx] = ARIFOMA_apply_all_schemes( ...
        corner, radarCtx, videok, numFrame, tIndex)
% ARIFOMA_APPLY_ALL_SCHEMES
%   Apply the active mitigation scheme to the victim and all interferers
%   for one corner view at one CPI micro-step.
%
% PURPOSE
% -------
%   This function is the scheme-dispatch layer of the ARIFOMA per-corner
%   CPI pipeline.  It updates the carrier frequency F0 of the victim sensor
%   and all interferer sensors according to their individually assigned
%   mitigation schemes (CONVENTIONAL, RANDOM_FH_PER_CPI, BATSINSPIRED, etc.).
%
%   The function reads the last-known F0 from the radar context (radarCtx),
%   delegates scheme logic to the victim and interferer scheme handlers, and
%   writes the updated F0 back to both the corner struct and the context.
%
%   After this function returns, corner.V.F0 and each corner.I(j).F0
%   reflect the carrier frequencies that will be used for waveform
%   generation in the current CPI micro-step.
%
% PIPELINE POSITION
% -----------------
%   Called from runCPI_perCorner (main_sim.m) as Stage 3 of the per-corner
%   CPI processing pipeline:
%
%     Stage 2  →  ARIFOMA_apply_all_schemes  →  Stage 4
%   (corner select)      (this function)      (build engine inputs)
%
% INPUTS
%   corner    (struct)  — per-corner scene bundle produced by
%                         get_corner_inputs; contains:
%                           corner.V   — victim parameter row
%                           corner.I   — interferer parameter rows (array)
%                           corner.ctx — current corner context
%   radarCtx  (struct)  — accumulated radar context struct; must contain:
%                           radarCtx.sensorSeedLUT — per-sensor seed / F0
%                           radarCtx.lastF0BySensor — last carrier frequency
%                                                     indexed by sensor key
%   videok    (int)     — current CPI micro-step index within the acquisition
%                         snapshot (1-based); passed to scheme handlers for
%                         intra-CPI adaptation logic
%   numFrame  (int)     — current frame index within the CPI block (1-based)
%   tIndex    (int)     — current acquisition time index (1-based);
%                         used by scheme handlers that accumulate history
%                         across snapshots
%
% OUTPUTS
%   corner    (struct)  — updated corner bundle with fields added/modified:
%                           corner.V            — victim with updated F0
%                           corner.I            — interferers with updated F0
%                           corner.ctx          — updated corner context
%                           corner.fv_base      — victim base carrier freq [Hz]
%                           corner.fi_base      — interferer base carrier
%                                                 frequencies (numI×1) [Hz]
%                           corner.schemeInterferer — active interferer scheme
%                                                     string (first interferer)
%   radarCtx  (struct)  — updated context; radarCtx.sensorSeedLUT.(keyV).f0
%                         is refreshed with the victim's new F0
%
% NOTES
% -----
%   - Interferer scheme: the current implementation uses the scheme of the
%     first interferer (schemeInterferer = schemeInterferer(1)) for all
%     interferers in this corner view.  This matches the homogeneous
%     deployment assumption used in the published case study.  For
%     heterogeneous interferer schemes, ARIFOMA_apply_interferer_schemes
%     would need to be called per-interferer.
%   - When no interferers are present (numI = 0), fi_base is returned as
%     an empty array and schemeInterferer is set to "NOTHING".
%   - The victim's F0 seed is read from ctx.lastF0BySensor.(keyV) so that
%     frequency-hopping schemes have access to the previous CPI's carrier.
%
% SEE ALSO
%   ARIFOMA_apply_victim_scheme,
%   ARIFOMA_apply_interferer_schemes,
%   runCPI_perCorner (main_sim.m)

% =========================================================================

%% ── Unpack Corner Bundle ─────────────────────────────────────────────────
V   = corner.V;
I   = corner.I;
ctx = corner.ctx;

%% ── Victim Scheme Application ────────────────────────────────────────────
%   Read the victim's last-known carrier frequency from the context.
%   This provides the starting point for frequency-hopping schemes that
%   need to know the previous CPI's F0 before deciding the next hop.

keyV     = sidKey(V.SensorID);
fv_base  = double(ctx.lastF0BySensor.(keyV));

[ctx, V, fv_base] = ARIFOMA_apply_victim_scheme( ...
    ctx, V, fv_base, ...
    string(V.Scheme), string(V.WaveformClean), V.Dir, ...
    V.bandMin, V.bandMax, V.B, V.FreqIncrHz, V.Nch, ...
    videok, numFrame, tIndex);

% Persist the victim's updated F0 to the seed LUT for use in future CPIs.
radarCtx.sensorSeedLUT.(keyV).f0 = V.F0;

%% ── Interferer Scheme Application ────────────────────────────────────────
%   All interferers in this corner view are processed together.
%   The scheme of the first interferer is applied to all (homogeneous
%   deployment assumption — see Notes).

numI = numel(I);

if numI > 0

    fi_base          = vertcat(I.F0);
    schemeInterferer = string({I.Scheme}).';
    schemeInterferer = schemeInterferer(1);   % first interferer scheme for all

    [ctx, I, fi_base] = ARIFOMA_apply_interferer_schemes( ...
        ctx, I, fi_base, ...
        schemeInterferer, string({I.Waveform}).', vertcat(I.Dir), ...
        vertcat(I.bandMin), vertcat(I.bandMax), vertcat(I.B), ...
        vertcat(I.FreqIncrHz), vertcat(I.Nch), ...
        videok, numFrame, tIndex, radarCtx.sensorSeedLUT);

else
    % No interferers present in this corner view.
    fi_base          = [];
    schemeInterferer = "NOTHING";
end

%% ── Pack Updated Corner Bundle ───────────────────────────────────────────
%   Write all updated objects back into the corner struct so that
%   ARIFOMA_build_engine_inputs receives the current-CPI carrier frequencies.

corner.ctx              = ctx;
corner.V                = V;
corner.I                = I;
corner.fv_base          = fv_base;
corner.fi_base          = fi_base;
corner.schemeInterferer = schemeInterferer;

end   % ARIFOMA_apply_all_schemes

% =========================================================================
%% ── Local Helper Functions ───────────────────────────────────────────────
% =========================================================================

% -------------------------------------------------------------------------
function key = sidKey(sid)
% SIDKEY  Convert a numeric sensor ID to its struct field key string.
%
%   MATLAB struct field names cannot begin with a digit, so sensor IDs are
%   prefixed with 'S' to form valid identifiers.
%
%   INPUT
%     sid  (double) — numeric sensor ID
%
%   OUTPUT
%     key  (char) — field key string, e.g. sensor ID 1234 → 'S1234'

    key = sprintf('S%d', round(double(sid)));
end