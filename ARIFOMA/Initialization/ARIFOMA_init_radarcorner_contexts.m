function radarCtx = ARIFOMA_init_radarcorner_contexts(radarCtx,cornerList)
% ARIFOMA_INIT_RADARCORNER_CONTEXTS
% This function creates the per-corner persistent state container used by the
% simulation/avoidance loop (e.g., do_tile()). The output is a struct with one
% field per requested corner label:
%
%   radarCtx.FL, radarCtx.FR, radarCtx.RL, radarCtx.RR, ...
%
% Each radarCtx.(corner) contains:
%   Core persistence
%     - lastF0BySensor : last center/start frequency per sensor (for return/jump policies)
%     - randPerFrameV/I: cached random choices for RAND_PER_FRAME(/CPI) schemes
%     - segV           : segmentation permutation state (SEGMENTED/SEGMENTED_PERM)
%     - bats_*         : BATS-inspired state (r_hat, mode, epsilon)
%     - dataLog        : minimal logging container (grown during runtime)
%
% Inputs
%   cornerList : string/cell array of corner labels, e.g. ["FL","FR","RL","RR"]
%
% Output
%   radarCtx   : struct containing initialized context per corner

%% -------------------- Default avoidance memory knobs --------------------
% Grouped knobs used by the avoidance strategies. Keep these centralized so
% experiments remain reproducible.
    for k = 1:numel(cornerList)
        cn = char(cornerList(k));
        radarCtx.(cn) = struct( ...
            'pairs', struct(),...
            'frameHadInterfV', false, ...
            'state', struct('victim',struct('has',false),'I',struct([])), ...
            't0', [],...
            'dataLog',  struct([]), ... % <<< seed correct shape
            'lastF0BySensor', struct(), ...
            'bats_r_hat', NaN, ...
            'bats_last_mode_V', "init", ...
            'bats_epsilon', 0.07);
end

end

