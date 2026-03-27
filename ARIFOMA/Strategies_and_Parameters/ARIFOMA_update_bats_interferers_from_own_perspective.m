function [ctx, info] = ARIFOMA_update_bats_interferers_from_own_perspective(ctx, V, I, perspectives, cfg)

info = struct([]);

for jj = 1:numel(I)

    if ~strcmpi(string(I(jj).Scheme), "BATSINSPIRED_FH_PER_CPI")
        continue;
    end

    keyIj = sprintf('S%d', I(jj).SensorID);

    if ~isfield(ctx,'waveformBySensor') || ~isfield(ctx.waveformBySensor, keyIj)
        continue;
    end

    Wi = ctx.waveformBySensor.(keyIj);

    % Build others for THIS interferer as ego
    [~, others_t, others_f, others_idx] = build_others_for_sensor( ...
        ctx, V, I, jj, perspectives);

    [t_i_all, f_i_all, idx_i_all] = pack_others_for_detector( ...
        others_t, others_f, others_idx);

    [~, ~, ~, ~, ~, pint, ~, ~, r_hat, f_int_min, f_int_max] = ...
        ARIFOMA_analytic_detector_batsinspired( ...
            Wi.t, Wi.f, Wi.idx, ...
            t_i_all, f_i_all, idx_i_all, ...
            cfg.BADC_Hz, cfg.boolIgnoreFlybackInterferer);

    % Ensure state structs exist
    if ~isfield(ctx,'policyBySensor') || ~isstruct(ctx.policyBySensor)
        ctx.policyBySensor = struct();
    end
    if ~isfield(ctx.policyBySensor, keyIj) || ~isstruct(ctx.policyBySensor.(keyIj))
        ctx.policyBySensor.(keyIj) = struct();
    end
    if ~isfield(ctx,'lastF0BySensor') || ~isstruct(ctx.lastF0BySensor)
        ctx.lastF0BySensor = struct();
    end

    polIj = ctx.policyBySensor.(keyIj);

    if isfinite(r_hat) && isfinite(f_int_min) && isfinite(f_int_max)
        polIj.bats_r_hat = r_hat;
        polIj.bats_next_f0 = bats_choose_next_f0_paper( ...
            I(jj).B, ...
            I(jj).Dir, ...
            I(jj).bandMin, ...
            I(jj).bandMax, ...
            r_hat, ...
            cfg.bats_eps, ...
            f_int_min, ...
            f_int_max);

        ctx.lastF0BySensor.(keyIj) = polIj.bats_next_f0;
    end

    ctx.policyBySensor.(keyIj) = polIj;

    info(end+1).SensorID = I(jj).SensorID; %#ok<AGROW>
    info(end).r_hat = r_hat;
    info(end).f_int_min = f_int_min;
    info(end).f_int_max = f_int_max;
    info(end).interference_probability = pint;
end

end

function [t_i, f_i, idx_i] = pack_others_for_detector(others_t, others_f, others_idx)

t_i   = [];
f_i   = [];
idx_i = [];

idxOffset = 0;

for m = 1:numel(others_t)
    tm = others_t{m};
    fm = others_f{m};
    im = others_idx{m};

    if isempty(tm) || isempty(fm) || isempty(im)
        continue;
    end

    tm = tm(:);
    fm = fm(:);
    im = im(:) + idxOffset;

    t_i   = [t_i; tm]; %#ok<AGROW>
    f_i   = [f_i; fm]; %#ok<AGROW>
    idx_i = [idx_i; im]; %#ok<AGROW>

    idxOffset = max(idx_i) + 1000;
end

end