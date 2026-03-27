function [ctx, V, fv_base] = ARIFOMA_apply_victim_scheme(ctx, V, fv_base, ...
    schemeVictim, waveformCleanV, dir_V, ...
    bandMin_V, bandMax_V, Bv, freqIncr_V, Nch, ...
    videok, nFrame, tIndex)
%APPLY_VICTIM_SCHEME
% Updates victim start frequency V.F0 according to schemeVictim.
%
% Requires these ctx fields to already exist:
%   ctx.randPerFrameV.tIndexStamp, ctx.randPerFrameV.f0
%   ctx.frameStampV, ctx.wasInterferedPrevFrameV, ctx.frameHadInterfV
%   ctx.bats_r_hat (only for BATSINSPIRED_PER_FRAME)

useStepped = strcmpi(string(waveformCleanV), 'STEPPED');

switch upper(string(schemeVictim))

    case "RAND_PER_FRAME"
        newBlock = (videok == 1) || (ctx.randPerFrameV.tIndexStamp ~= tIndex);
        if newBlock
            ctx.randPerFrameV.tIndexStamp = tIndex;
            ctx.randPerFrameV.f0 = pickStart(bandMin_V, bandMax_V, Bv, freqIncr_V, Nch, dir_V, useStepped);
        end
        fv_base = ctx.randPerFrameV.f0;
        V.F0 = fv_base;

    case "RAND_PER_FRAME_INTERFERED"
        newFrame = (videok == 1) || (ctx.frameStampV ~= tIndex);
        if newFrame
            if ctx.wasInterferedPrevFrameV
                fv_base = pickStart(bandMin_V, bandMax_V, Bv, freqIncr_V, Nch, dir_V, useStepped);
            end
            ctx.frameStampV = tIndex;
            ctx.frameHadInterfV = false;
        end
        V.F0 = fv_base;

    case "BATSINSPIRED_PER_FRAME"
        cpiLen = max(1, floor(nFrame/2));
        isCpiStart = (videok == 1) || (videok == cpiLen + 1);
        if isCpiStart
            r_hat = double(ctx.bats_r_hat);  % must exist; may be NaN (your choice)
            epsBat = double(ctx.bats_epsilon);

            f0_new = fv_base;

            if isfinite(r_hat)
                f_intc = fv_base + r_hat*Bv;

                if abs(r_hat - 0.5) < epsBat
                    f0_new = pickStart(bandMin_V, bandMax_V, Bv, freqIncr_V, Nch, dir_V, useStepped);
                    ctx.bats_last_mode_V = "random";
                elseif r_hat < 0.5
                    if dir_V >= 0
                        f0_new = min(f_intc, bandMax_V - Bv);
                    else
                        f0_new = min(bandMax_V, f_intc + Bv/2);
                    end
                    ctx.bats_last_mode_V = "shift_up";
                else
                    if dir_V >= 0
                        f0_new = max(bandMin_V, f_intc - Bv);
                    else
                        f0_new = max(bandMin_V + Bv, f_intc);
                    end
                    ctx.bats_last_mode_V = "shift_down";
                end
            else
                ctx.bats_last_mode_V = "init";
            end

            fv_base = f0_new;
            V.F0 = f0_new;
        else
            V.F0 = fv_base;
        end

    otherwise
        % default: keep fv_base
        V.F0 = fv_base;
end

end
