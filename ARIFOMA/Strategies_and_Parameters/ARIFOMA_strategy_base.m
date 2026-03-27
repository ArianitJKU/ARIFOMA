function out = ARIFOMA_strategy_base(engine)

c = engine.corner;
V = c.V;
I = c.I;

outS = ARIFOMA_conventional_scene_runCPI( ...
    engine.ax, engine.numCPI, engine.cfg, c.cn, engine.ctx, V, I, ...
    c.fv_base, c.fi_base, c.schemeInterferer);

% update adaptive interferers from THEIR perspective
[ctx, batsInfo] = ARIFOMA_update_bats_interferers_from_own_perspective( ...
    outS.ctx, V, I, c.perspectives, engine.cfg);

out.ctx = ctx;
out.metrics = outS;
out.metrics.batsInterferers = batsInfo;
out.score_trace = struct();
end
