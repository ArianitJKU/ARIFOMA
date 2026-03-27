function out = ARIFOMA_run_strategy(engine)

scheme = upper(string(engine.corner.V.Scheme));

switch scheme
    case {"CONVENTIONAL","FAR","RANDOM_FH_PER_CPI","SEGMENTED"}
        out = ARIFOMA_strategy_base(engine);

    case "BATSINSPIRED_FH_PER_CPI"
        out = ARIFOMA_strategy_batsinspiredfhperCPI(engine);

    otherwise
        error("Unknown victim scheme/strategy: %s", string(engine.corner.V.Scheme));
end
