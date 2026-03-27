function ctx = ARIFOMA_update_interferer_bats_stats_end_cpi(ctx, schemeInterferer, I, videok, nFrame)
%UPDATE_INTERFERER_BATS_STATS_END_CPI
% Finalize per-interferer CPI-mean r_hat from memBySensor accumulators.

if ~( strcmpi(string(schemeInterferer),'BATSINSPIRED_PER_FRAME') || ...
      strcmpi(string(schemeInterferer),'RAND_PER_FRAME_INTERFERED') )
    return;
end

cpiLen = max(1, floor(nFrame/2));
isEnd1 = (videok == cpiLen);
isEnd2 = (videok == 2*cpiLen);
if ~(isEnd1 || isEnd2)
    return;
end

nI = numel(I);
for j = 1:nI
    keyI = sidKey(I(j).SensorID);
    memI = ctx.memBySensor.(keyI);

    % Your accumulators must exist and be valid numbers.
    memI.bats_r_hat = memI.bats_sum_r_hat / memI.bats_n_subframes;

    % Reset accumulators for next CPI
    memI.bats_sum_r_hat = 0;
    memI.bats_n_subframes = 0;

    ctx.memBySensor.(keyI) = memI;
end

end
