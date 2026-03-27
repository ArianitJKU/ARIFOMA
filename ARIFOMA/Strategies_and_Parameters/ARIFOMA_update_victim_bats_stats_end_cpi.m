function ctx = ARIFOMA_update_victim_bats_stats_end_cpi(ctx, schemeVictim, tIndex, videok, nFrame)
%UPDATE_VICTIM_BATS_STATS_END_CPI
% At CPI end, compute mean of frame(k).centerofIntf over that CPI and store in ctx.bats_r_hat.

if ~strcmpi(string(schemeVictim), 'BATSINSPIRED_PER_FRAME')
    return;
end

cpiLen = max(1, floor(nFrame/2));
isEnd1 = (videok == cpiLen);
isEnd2 = (videok == 2*cpiLen);

if ~(isEnd1 || isEnd2)
    return;
end

if isEnd1
    k0 = 1;       k1 = cpiLen;
else
    k0 = cpiLen+1; k1 = 2*cpiLen;
end

rvals = [];
for k = k0:k1
    Fr = ctx.log.step(tIndex).frame(k);
    if isfield(Fr,'centerofIntf') && isfinite(Fr.centerofIntf)
        rvals(end+1) = double(Fr.centerofIntf); %#ok<AGROW>
    end
end

if ~isempty(rvals)
    ctx.bats_r_hat = mean(rvals);
end

end
