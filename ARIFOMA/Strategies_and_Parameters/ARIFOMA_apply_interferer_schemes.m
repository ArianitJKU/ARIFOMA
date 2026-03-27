function [ctx, I, fi_base] = ARIFOMA_apply_interferer_schemes(ctx, I, fi_base, ...
        schemeInterferer, waveformCleanI, dir_I, ...
        bandMin_I, bandMax_I, Bi, freqIncr_I, NchI, ...
        videok, nFrame, tIndex, sensorSeedLUT)
% ARIFOMA_APPLY_INTERFERER_SCHEMES
%   Update carrier frequencies of all interferer sensors for one CPI micro-step.
%
% PURPOSE
% -------
%   Applies the active mitigation scheme to all interferers in a corner view,
%   updating I(j).F0 and fi_base(j) for each interferer j according to the
%   scheme assigned in schemeInterferer.
%
%   Three scheme families are supported:
%     RAND_PER_FRAME            — random F0 drawn once per acquisition block,
%                                 shared across all interferers in the view
%     RAND_PER_FRAME_INTERFERED — random F0 drawn per-interferer at CPI
%                                 boundaries, conditioned on interference
%                                 history in ctx.memBySensor
%     BATSINSPIRED_PER_FRAME    — BATS-inspired directional frequency shift
%                                 derived from the interference centroid
%                                 estimate r_hat stored in ctx.memBySensor
%     (otherwise)               — no update; F0 remains as initialised
%
% CPI BOUNDARY DETECTION
% ----------------------
%   Schemes that update only at CPI boundaries use isCpiStart:
%     isCpiStart = true  when videok == 1  (first micro-step of snapshot)
%                     OR when videok == cpiLen + 1 (start of second CPI half)
%   where cpiLen = floor(nFrame / 2).
%
% INPUTS
%   ctx             (struct)       — current corner context; required fields
%                                   depend on active scheme (see Notes)
%   I               (struct array) — interferer parameter structs; each must
%                                    contain: SensorID, F0, Scheme, Waveform
%   fi_base         (double, nI×1) — current interferer base carrier
%                                    frequencies [Hz]; updated in-place
%   schemeInterferer (string)      — scheme identifier applied to all
%                                    interferers in this view (from first
%                                    interferer; see ARIFOMA_apply_all_schemes)
%   waveformCleanI  (string array, nI×1) — canonical waveform type per
%                                    interferer ("FMCW" or "STEPPED FMCW")
%   dir_I           (double, nI×1) — chirp slope sign per interferer
%                                    (+1 positive, -1 negative)
%   bandMin_I       (double, nI×1) — RF band lower bound per interferer [Hz]
%   bandMax_I       (double, nI×1) — RF band upper bound per interferer [Hz]
%   Bi              (double, nI×1) — sweep bandwidth per interferer [Hz]
%   freqIncr_I      (double, nI×1) — stepped frequency increment per
%                                    interferer [Hz]; 0 for standard FMCW
%   NchI            (double, nI×1) — number of chirps per CPI per interferer
%   videok          (int)          — current CPI micro-step index (1-based)
%   nFrame          (int)          — total CPI micro-steps per acquisition
%   tIndex          (int)          — current acquisition time index (1-based)
%   sensorSeedLUT   (struct)       — per-sensor seed / F0 lookup table;
%                                    indexed by sensor key S<id>
%
% OUTPUTS
%   ctx      (struct)       — updated context; modified fields depend on
%                             active scheme:
%                               ctx.lastF0BySensor.(keyI)
%                               ctx.memBySensor.(keyI)
%                               ctx.state.I(j).f0
%                               ctx.randPerFrameI  (RAND_PER_FRAME only)
%   I        (struct array) — interferer structs with I(j).F0 updated
%   fi_base  (double, nI×1) — updated interferer base carrier frequencies [Hz]
%
% NOTES
% -----
%   Context fields required per scheme:
%
%   RAND_PER_FRAME:
%     ctx.randPerFrameI.tIndexStamp — acquisition index of last F0 draw
%     ctx.randPerFrameI.f0          — last randomly drawn F0
%
%   RAND_PER_FRAME_INTERFERED:
%     ctx.memBySensor.(keyI).bats_r_hat — interference centroid estimate;
%                                         if finite, triggers a new F0 draw
%
%   BATSINSPIRED_PER_FRAME:
%     ctx.memBySensor.(keyI).bats_r_hat — interference centroid estimate r_hat
%     ctx.bats_epsilon                  — convergence threshold ε;
%                                         |r_hat - 0.5| < ε triggers random hop
%
%   RAND_PER_FRAME applies the same base F0 to all interferers, with
%   per-interferer clamping to ensure F0 stays within the RF band and does
%   not cause the sweep to exceed bandMax_I.
%
%   The sensorSeedLUT write inside RAND_PER_FRAME is a no-op in the current
%   architecture because MATLAB passes structs by value — the caller's LUT
%   is not modified.  This is noted in-code and flagged for future refactor.
%
% SEE ALSO
%   ARIFOMA_apply_all_schemes,
%   ARIFOMA_apply_victim_scheme,
%   ARIFOMA_analytic_detector_batsinspired

% =========================================================================

nI = numel(I);
if nI == 0
    return;
end

%% ── CPI Boundary Detection ───────────────────────────────────────────────
%   isCpiStart is true at the first micro-step of each CPI half.
%   Schemes that update only at CPI boundaries gate on this flag.

cpiLen     = max(1, floor(nFrame / 2));
isCpiStart = (videok == 1) || (videok == cpiLen + 1);

%% ── Scheme Dispatch ──────────────────────────────────────────────────────

switch upper(string(schemeInterferer))

    % ── RAND_PER_FRAME ────────────────────────────────────────────────────
    %   Draw one random F0 per acquisition block (shared across all
    %   interferers in this view).  A new draw is triggered on the first
    %   micro-step of the snapshot or when the acquisition index changes.
    %   All interferers receive the same base F0, individually clamped to
    %   their RF band boundaries.

    case "RAND_PER_FRAME"

        % Determine whether this is the start of a new acquisition block.
        newBlock = (videok == 1) || ...
                   (ctx.randPerFrameI.tIndexStamp ~= tIndex);

        if newBlock
            useStepped0 = strcmpi(string(waveformCleanI(1)), 'STEPPED');
            ctx.randPerFrameI.tIndexStamp = tIndex;
            ctx.randPerFrameI.f0 = pickStart( ...
                bandMin_I(1), bandMax_I(1), Bi(1), ...
                freqIncr_I(1), NchI(1), dir_I(1), useStepped0);
        end

        baseF0 = double(ctx.randPerFrameI.f0);

        % Apply the shared base F0 to each interferer with individual
        % band / sweep-span clamping.
        for j = 1:nI
            f0j = baseF0;

            if dir_I(j) >= 0
                % Positive slope: F0 is the sweep start; clamp so sweep
                % does not exceed bandMax.
                f0j = min(max(f0j, bandMin_I(j)), bandMax_I(j) - Bi(j));
            else
                % Negative slope: F0 is the sweep end; clamp accordingly.
                f0j = max(min(f0j, bandMax_I(j)), bandMin_I(j) + Bi(j));
            end

            fi_base(j) = f0j;
            I(j).F0    = f0j;

            % NOTE: sensorSeedLUT is passed by value in MATLAB — this write
            % does not propagate to the caller.  Retained for future refactor
            % to a handle-based LUT.
            keyI = sidKey(I(j).SensorID);
            sensorSeedLUT.(keyI).f0 = f0j;   %#ok<NASGU>
        end

    % ── RAND_PER_FRAME_INTERFERED ─────────────────────────────────────────
    %   Draw a new random F0 per-interferer at CPI boundaries, but only
    %   when the interferer's interference history indicates it has been
    %   detected (bats_r_hat is finite in ctx.memBySensor).
    %   If no history exists, F0 is left unchanged.

    case "RAND_PER_FRAME_INTERFERED"

        if isCpiStart
            for j = 1:nI
                keyI  = sidKey(I(j).SensorID);
                memI  = ctx.memBySensor.(keyI);

                useStepped = strcmpi(string(waveformCleanI(j)), 'STEPPED');

                f0_new = double(I(j).F0);

                % Redraw only if interference has been previously detected.
                if isfield(memI, 'bats_r_hat') && isfinite(memI.bats_r_hat)
                    f0_new = pickStart( ...
                        bandMin_I(j), bandMax_I(j), Bi(j), ...
                        freqIncr_I(j), NchI(j), dir_I(j), useStepped);
                end

                I(j).F0    = f0_new;
                fi_base(j) = f0_new;

                ctx.lastF0BySensor.(keyI) = f0_new;
                ctx.state.I(j).F0         = f0_new;
            end
        end

    % ── BATSINSPIRED_PER_FRAME ────────────────────────────────────────────
    %   BATS-inspired directional frequency shift applied at CPI boundaries.
    %   Uses the interference centroid estimate r_hat from ctx.memBySensor
    %   to decide whether to shift F0 up, shift F0 down, or hop randomly:
    %
    %     |r_hat - 0.5| < ε  →  random hop  (interference centred; uncertain)
    %     r_hat < 0.5        →  shift up     (interference in lower chirp half)
    %     r_hat >= 0.5       →  shift down   (interference in upper chirp half)
    %
    %   If no r_hat is available (first CPI or no prior detection), F0 is
    %   left unchanged and the mode is recorded as "init".

    case "BATSINSPIRED_PER_FRAME"

        if isCpiStart
            epsBat = double(ctx.bats_epsilon);

            for j = 1:nI
                keyI  = sidKey(I(j).SensorID);
                memI  = ctx.memBySensor.(keyI);

                useStepped = strcmpi(string(waveformCleanI(j)), 'STEPPED');

                f0_new = double(I(j).F0);

                if isfield(memI, 'bats_r_hat') && isfinite(memI.bats_r_hat)

                    r_hat  = double(memI.bats_r_hat);
                    f_intc = f0_new + r_hat * double(Bi(j));

                    if abs(r_hat - 0.5) < epsBat
                        % ── Random hop: centroid near midpoint → uncertain ──
                        f0_new = pickStart( ...
                            bandMin_I(j), bandMax_I(j), Bi(j), ...
                            freqIncr_I(j), NchI(j), dir_I(j), useStepped);
                        memI.bats_last_mode_I = "random";

                    elseif r_hat < 0.5
                        % ── Shift up: interference in lower chirp half ─────
                        if dir_I(j) >= 0
                            f0_new = min(f_intc, bandMax_I(j) - Bi(j));
                        else
                            f0_new = min(bandMax_I(j), f_intc + Bi(j)/2);
                        end
                        memI.bats_last_mode_I = "shift_up";

                    else
                        % ── Shift down: interference in upper chirp half ───
                        if dir_I(j) >= 0
                            f0_new = max(bandMin_I(j), f_intc - Bi(j));
                        else
                            f0_new = max(bandMin_I(j) + Bi(j), f_intc);
                        end
                        memI.bats_last_mode_I = "shift_down";
                    end

                else
                    % No prior r_hat available — leave F0 unchanged.
                    memI.bats_last_mode_I = "init";
                end

                I(j).F0    = f0_new;
                fi_base(j) = f0_new;

                ctx.lastF0BySensor.(keyI) = f0_new;
                ctx.memBySensor.(keyI)    = memI;
                ctx.state.I(j).F0         = f0_new;
            end
        end

    otherwise
        % ── No update: keep fi_base and I(j).F0 as initialised ───────────
        %   Covers "CONVENTIONAL" and any unrecognised scheme strings.

end   % scheme switch

end   % ARIFOMA_apply_interferer_schemes

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