function out = ARIFOMA_strategy_batsinspiredfhperCPI(engine)

c   = engine.corner;
V   = c.V;
I   = c.I;
ctx = engine.ctx;

keyV = sprintf('S%d', V.SensorID);

% Make sure state structs exist
if ~isfield(ctx,'policyBySensor') || ~isstruct(ctx.policyBySensor)
    ctx.policyBySensor = struct();
end
if ~isfield(ctx.policyBySensor, keyV) || ~isstruct(ctx.policyBySensor.(keyV))
    ctx.policyBySensor.(keyV) = struct();
end
if ~isfield(ctx,'lastF0BySensor') || ~isstruct(ctx.lastF0BySensor)
    ctx.lastF0BySensor = struct();
end

polV = ctx.policyBySensor.(keyV);

% Use the F0 scheduled from the PREVIOUS CPI for THIS sensor
fv_this = c.fv_base;
if isfield(polV,'bats_next_f0') && isfinite(polV.bats_next_f0)
    fv_this = polV.bats_next_f0;
elseif isfield(ctx.lastF0BySensor, keyV) && isfinite(ctx.lastF0BySensor.(keyV))
    fv_this = ctx.lastF0BySensor.(keyV);
end

% Run this CPI from the victim's own perspective
outS = ARIFOMA_batsinspired_scene_runCPI( ...
    engine.ax, engine.numCPI,engine.cfg, c.cn, ctx, V, I, ...
    fv_this, c.fi_base, c.schemeInterferer);

ctx = outS.ctx;

% --------------------------------------------------
% (1) Update bats state for THIS victim sensor
% --------------------------------------------------
if isfield(outS,'all_r_hat_sub') && isfinite(outS.all_r_hat_sub) && ...
   isfield(outS,'f_int_min')     && isfinite(outS.f_int_min)     && ...
   isfield(outS,'f_int_max')     && isfinite(outS.f_int_max)

    polV.bats_r_hat = outS.all_r_hat_sub;

    polV.bats_next_f0 = bats_choose_next_f0_paper( ...
        V.B, ...
        V.Dir, ...
        V.bandMin, ...
        V.bandMax, ...
        polV.bats_r_hat, ...
        engine.cfg.bats_eps, ...
        outS.f_int_min, ...
        outS.f_int_max);

else
    % no valid estimate -> keep previous frequency
    polV.bats_next_f0 = fv_this;
end

ctx.policyBySensor.(keyV) = polV;
ctx.lastF0BySensor.(keyV) = polV.bats_next_f0;

% --------------------------------------------------
% (2) ALSO update any interferers that are bats-inspired
%     from THEIR own perspective
% --------------------------------------------------
[ctx, batsInfo] = ARIFOMA_update_bats_interferers_from_own_perspective( ...
    ctx, V, I, c.perspectives, engine.cfg);

outS.ctx = ctx;

out.ctx = ctx;
out.metrics = outS;
out.metrics.batsVictim = struct( ...
    'SensorID', V.SensorID, ...
    'r_hat', getfield_safe(polV,'bats_r_hat',NaN), ...
    'next_f0', getfield_safe(polV,'bats_next_f0',NaN));
out.metrics.batsInterferers = batsInfo;
out.score_trace = struct();
end
function f0_next = bats_choose_next_f0_paper( ...
    B, dir, bandMin, bandMax, rhat, eps, f_int_min, f_int_max)

    % Random shift if no clear decision can be made
    if abs(rhat - 0.5) < eps
        f0_next = randomF0ForNormal(B, dir, bandMin, bandMax);
        return;
    end

    % Paper rule:
    % rhat < 0.5  -> upward shift -> place above highest interfered frequency
    % rhat > 0.5  -> downward shift -> place below lowest interfered frequency
    if rhat < 0.5
        if dir >= 0
            cand = f_int_max;
        else
            cand = f_int_max + B;
        end
    else
        if dir >= 0
            cand = f_int_min - B;
        else
            cand = f_int_min;
        end
    end

    f0_next = clampF0(cand, B, dir, bandMin, bandMax);
end
function val = getfield_safe(S, fld, defaultVal)
    if isstruct(S) && isfield(S,fld) && ~isempty(S.(fld)) && isfinite(S.(fld))
        val = S.(fld);
    else
        val = defaultVal;
    end
end
function f0r = randomF0ForNormal(B, dir, bandMin, bandMax)

    if dir >= 0
        lo = bandMin;
        hi = bandMax - B;
    else
        lo = bandMin + B;
        hi = bandMax;
    end

    if hi < lo
        error('randomF0ForNormal: no valid f0 range inside band.');
    end

    if hi == lo
        f0r = lo;
    else
        f0r = lo + rand*(hi - lo);
    end
end
function f0c = clampF0(f0, B, dir, bandMin, bandMax)

    if dir >= 0
        f0c = min(max(f0, bandMin), bandMax - B);
    else
        f0c = min(max(f0, bandMin + B), bandMax);
    end
end