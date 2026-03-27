function [dt_all, df_all, t_start_all, t_end_all, mask_chirp, ...
          interference_probability, frac_chirp, intervals_local_per_chirp, ...
          r_hat, f_int_min, f_int_max] = ...
    ARIFOMA_analytic_detector_batsinspired(t_v, f_v, idx_v, ...
                                           t_i, f_i, idx_i, ...
                                           BADC, boolIgnoreFlybackInterferer)
% ARIFOMA_ANALYTIC_DETECTOR_BATSINSPIRED
%   Analytic FMCW interference interval detector with BATS-inspired observables.
%
% PURPOSE
% -------
%   Detects all time intervals during which an interfering FMCW radar falls
%   within the victim's low-pass filter (LPF) bandwidth, and computes the
%   unified interference-related figures of merit (FOMs) defined in ARIFOMA.
%
%   The detector is fully analytic: each chirp is modelled as a set of
%   piecewise-linear monotonic segments and the condition
%   |f_V(t) - f_I(t)| <= BADC is solved in closed form.  This avoids all
%   sampling-grid artifacts that affect time-domain detectors.
%
%   In addition to the standard ARIFOMA FOMs, this function computes the
%   BATS-inspired observables r_hat, f_int_min, and f_int_max used by the
%   BATSINSPIRED_FH_PER_CPI mitigation scheme to estimate the centre of
%   interfered samples and redirect the victim's carrier frequency.
%
% KEY DESIGN DECISIONS
% --------------------
%   - Victim TX-only assumption: the victim is assumed to transmit only
%     during its main ramp (boolIgnoreFlybackVictim = true, hardcoded).
%     This reflects typical automotive FMCW practice.
%   - Interferer flyback: controlled by boolIgnoreFlybackInterferer.
%     Set true to model "TX only" interferers (main ramp only).
%     Set false to include flyback and return segments.
%   - Stepped FMCW: fully supported.  Frequency steps occur between chirps
%     (i.e. between segments), not inside the ramp, so piecewise-linear
%     segment modelling handles them correctly without modification.
%   - Interval merging: overlapping intervals from different interferer
%     segments are merged before metric computation to avoid double-counting.
%
% ALGORITHM OVERVIEW
% ------------------
%   1. Build victim chirp models: convert each victim chirp's (t,f) samples
%      into a single primary linear segment (main TX ramp).
%   2. Build interferer segment list: convert all interferer chirps into
%      monotonic linear segments (primary or all, per boolIgnoreFlybackInterferer).
%   3. For each victim chirp:
%        a. For each interferer segment, compute the time overlap window.
%        b. Model the IF frequency difference as D(t) = p*t + q (linear).
%        c. Solve |D(t)| <= BADC analytically over the overlap window.
%        d. Merge all resulting intervals.
%        e. Compute per-chirp FOMs and BATS observables.
%   4. Aggregate across all victim chirps.
%
% INPUTS
%   t_v   (double array, N×1) — victim instantaneous frequency sample times [s]
%   f_v   (double array, N×1) — victim instantaneous frequencies [Hz]
%   idx_v (double array, N×1) — victim chirp indices; NaN entries are ignored
%   t_i   (double array, M×1) — interferer instantaneous frequency sample times [s]
%   f_i   (double array, M×1) — interferer instantaneous frequencies [Hz]
%   idx_i (double array, M×1) — interferer chirp indices; NaN entries are ignored
%   BADC  (double)            — LPF / ADC bandwidth threshold [Hz];
%                               interference is detected when
%                               |f_V(t) - f_I(t)| <= BADC
%   boolIgnoreFlybackInterferer (logical)
%                             — true:  model interferer as TX-only (main ramp only)
%                               false: include flyback and return segments
%
% OUTPUTS
%   dt_all        (double array, 1×P) — duration of each detected overlap interval [s]
%   df_all        (double array, 1×P) — peak |f_V - f_I| within each interval [Hz];
%                                       conservatively set to BADC
%   t_start_all   (double array, 1×P) — absolute start time of each interval [s]
%   t_end_all     (double array, 1×P) — absolute end time of each interval [s]
%   mask_chirp    (logical array, 1×K)— true if victim chirp k has any overlap
%   interference_probability (double) — fraction of victim chirps with any overlap
%                                       (per-CPI interference ratio FOM)
%   frac_chirp    (double array, 1×K) — fraction of victim TX time overlapped
%                                       per chirp (normalised interference duration FOM)
%   intervals_local_per_chirp (cell, K×1)
%                             — per-chirp overlap intervals [tStart, tEnd] relative
%                               to the victim TX start time [s]
%   r_hat         (double)    — BATS interference centroid estimate;
%                               duration-weighted mean normalised position of
%                               interfered samples within the victim chirp [0, 1];
%                               NaN if no interference detected
%   f_int_min     (double)    — minimum victim RF frequency that is interfered [Hz];
%                               NaN if no interference detected
%   f_int_max     (double)    — maximum victim RF frequency that is interfered [Hz];
%                               NaN if no interference detected
%
% NOTES
% -----
%   - df_all is set conservatively to BADC for all intervals.  To obtain
%     exact peak frequency differences, the interferer segment responsible
%     for each interval would need to be stored during Step 3b — this is
%     a straightforward extension if needed.
%   - The function returns empty arrays (not errors) when no interference
%     is detected or when t_i is empty.
%   - r_hat is the continuous-time analogue of the discrete estimator in
%     the BATS paper [Bechter et al., IEEE MTT-S ICMIM 2016]:
%       r_hat = sum(u_bar * dur) / sum(dur)
%     where u_bar is the normalised midpoint of each interfered interval.
%
% REFERENCES
%   [1] J. Bechter, C. Sippel, and C. Waldschmidt, "Bats-inspired frequency
%       hopping for mitigation of interference between automotive radars,"
%       IEEE MTT-S ICMIM, 2016.
%
% SEE ALSO
%   ARIFOMA_run_strategy, ARIFOMA_apply_all_schemes,
%   ARIFOMA_build_engine_inputs

% =========================================================================

%% ── Input Sanitisation ───────────────────────────────────────────────────
%   Enforce column vectors for all sample arrays.

t_v   = t_v(:);    f_v   = f_v(:);    idx_v = idx_v(:);
t_i   = t_i(:);    f_i   = f_i(:);    idx_i = idx_i(:);
BADC  = double(BADC);

%% ── BATS Observable Accumulators ─────────────────────────────────────────
%   r_hat is computed as a duration-weighted mean of the normalised
%   interference centroid positions across all victim chirps.

r_num     = 0;       % Numerator:   sum(u_bar * dur)
r_den     = 0;       % Denominator: sum(dur)
r_hat     = NaN;     % Final estimate (NaN until at least one interval found)
f_int_min = NaN;     % Minimum interfered victim RF frequency [Hz]
f_int_max = NaN;     % Maximum interfered victim RF frequency [Hz]

%% ── Victim TX-Only Assumption ────────────────────────────────────────────
%   Automotive FMCW victims are assumed to transmit only during the main
%   ramp, not during the flyback or return interval.  This flag is
%   hardcoded to true and not exposed as a parameter to keep the interface
%   consistent with the ARIFOMA interference model.

boolIgnoreFlybackVictim = true;

% =========================================================================
%% ── Step 1: Build Victim Chirp Models ────────────────────────────────────
% =========================================================================
%   Each victim chirp is converted to a single primary linear segment
%   representing the main TX ramp.  The segment is stored as the linear
%   model f_V(t) = a*t + b over the interval [t0, t1].

chirp_ids_v = unique(idx_v(isfinite(idx_v)))';
K           = numel(chirp_ids_v);

victimSeg = repmat(struct('t0',NaN, 't1',NaN, 'a',NaN, 'b',NaN), K, 1);

for k = 1:K
    mk = isfinite(idx_v) & (idx_v == chirp_ids_v(k));
    tv = t_v(mk);
    fv = f_v(mk);

    if boolIgnoreFlybackVictim
        segs = chirp_to_segments(tv, fv, 'primary');   % main TX ramp only
    else
        segs = chirp_to_segments(tv, fv, 'all');        % all segments
    end

    if isempty(segs)
        continue;
    end

    % When 'all' mode is used, the first segment is still taken as the
    % primary TX ramp for chirp-duration normalisation.
    victimSeg(k) = segs(1);
end

% =========================================================================
%% ── Step 2: Build Interferer Segment List ────────────────────────────────
% =========================================================================
%   All interferer chirps are pre-converted to monotonic linear segments.
%   Pre-computing this list avoids redundant segment extraction inside the
%   per-victim-chirp loop.

interfererSegments = struct('t0',{}, 't1',{}, 'a',{}, 'b',{});

if ~isempty(t_i)
    chirp_ids_i = unique(idx_i(isfinite(idx_i)))';

    for j = 1:numel(chirp_ids_i)
        mj = isfinite(idx_i) & (idx_i == chirp_ids_i(j));
        tj = t_i(mj);
        fj = f_i(mj);

        if boolIgnoreFlybackInterferer
            segs = chirp_to_segments(tj, fj, 'primary');   % TX only
        else
            segs = chirp_to_segments(tj, fj, 'all');        % include flyback
        end

        if ~isempty(segs)
            interfererSegments = [interfererSegments; segs(:)];   %#ok<AGROW>
        end
    end
end

% =========================================================================
%% ── Step 3: Per-Victim-Chirp Analytic Interference Detection ─────────────
% =========================================================================
%   For each victim chirp, solve |f_V(t) - f_I(t)| <= BADC analytically
%   against every interferer segment whose time window overlaps the victim
%   TX interval.

dt_all      = [];
df_all      = [];
t_start_all = [];
t_end_all   = [];

mask_chirp                = false(1, K);
frac_chirp                = zeros(1, K);
intervals_local_per_chirp = cell(K, 1);

for k = 1:K

    v = victimSeg(k);

    % Skip chirps with undefined or degenerate TX intervals.
    if ~isfinite(v.t0) || ~isfinite(v.t1) || v.t1 <= v.t0
        intervals_local_per_chirp{k} = zeros(0, 2);
        continue;
    end

    t0V = v.t0;
    t1V = v.t1;
    Tv  = max(1e-12, t1V - t0V);   % victim TX duration; guard against zero

    % ── 3a: Collect candidate overlap intervals ───────────────────────────
    %   Test every interferer segment against this victim chirp.

    I_abs = zeros(0, 2);

    for s = 1:numel(interfererSegments)
        Iseg = interfererSegments(s);

        % Time window common to both the victim TX and this interferer segment.
        tL = max(t0V, Iseg.t0);
        tR = min(t1V, Iseg.t1);

        if ~(isfinite(tL) && isfinite(tR)) || (tR <= tL)
            continue;
        end

        % ── 3b: Analytic frequency-difference model ───────────────────────
        %   f_V(t) = v.a*t + v.b
        %   f_I(t) = Iseg.a*t + Iseg.b
        %   D(t)   = f_V(t) - f_I(t) = p*t + q  (linear in t)
        p = v.a - Iseg.a;
        q = v.b - Iseg.b;

        % ── 3c: Solve |p*t + q| <= BADC over [tL, tR] ────────────────────
        J = abs_linear_interval(p, q, BADC, tL, tR);   % returns [] or [ts te]

        if ~isempty(J)
            I_abs(end+1, :) = J;   %#ok<AGROW>
        end
    end

    if isempty(I_abs)
        intervals_local_per_chirp{k} = zeros(0, 2);
        continue;
    end

    % ── 3d: Merge overlapping intervals ───────────────────────────────────
    %   Multiple interferer segments may produce overlapping detections.
    %   Merging prevents double-counting in the FOM calculations.
    I_abs = merge_intervals(I_abs, 0);

    % ── 3e: BATS Observable Accumulation ──────────────────────────────────
    %   For each merged interval, accumulate the duration-weighted
    %   normalised midpoint u_bar and the interfered RF frequency range.

    for r = 1:size(I_abs, 1)
        ts  = I_abs(r, 1);
        te  = I_abs(r, 2);
        dur = te - ts;

        if dur <= 0
            continue;
        end

        % Normalised start and end positions within victim chirp [0, 1].
        ua = min(max((ts - t0V) / Tv, 0), 1);
        ub = min(max((te - t0V) / Tv, 0), 1);

        % Duration-weighted centroid accumulation.
        % For uniformly sampled ramps this matches the discrete BATS estimator.
        ubar  = 0.5 * (ua + ub);
        r_num = r_num + ubar * dur;
        r_den = r_den + dur;

        % Map interfered time window to victim RF frequencies using the
        % victim's linear ramp model: f = a*t + b.
        fa   = v.a * ts + v.b;
        fb   = v.a * te + v.b;
        f_lo = min(fa, fb);
        f_hi = max(fa, fb);

        if isnan(f_int_min)
            f_int_min = f_lo;
            f_int_max = f_hi;
        else
            f_int_min = min(f_int_min, f_lo);
            f_int_max = max(f_int_max, f_hi);
        end
    end   % BATS accumulation loop

    % ── Store local intervals (relative to victim TX start) ───────────────
    I_loc        = I_abs;
    I_loc(:, 1)  = I_abs(:, 1) - t0V;
    I_loc(:, 2)  = I_abs(:, 2) - t0V;
    intervals_local_per_chirp{k} = I_loc;

    % ── Per-chirp FOM computation ─────────────────────────────────────────
    mask_chirp(k) = true;

    % Normalised interference duration: total overlapped time / chirp TX time.
    dur_total  = sum(max(0, I_abs(:,2) - I_abs(:,1)));
    frac_chirp(k) = dur_total / Tv;

    % ── Append to global interval lists ───────────────────────────────────
    %   One entry per merged interval; df_all is conservatively set to BADC.
    %   (See NOTES for exact df extension.)
    for r = 1:size(I_abs, 1)
        ts = I_abs(r, 1);
        te = I_abs(r, 2);

        dt_all(end+1)      = max(0, te - ts);   %#ok<AGROW>
        df_all(end+1)      = BADC;               %#ok<AGROW>  conservative bound
        t_start_all(end+1) = ts;                 %#ok<AGROW>
        t_end_all(end+1)   = te;                 %#ok<AGROW>
    end

    % Update r_hat after each chirp so it reflects all chirps processed so far.
    if r_den > 0
        r_hat = r_num / r_den;
    end

end   % victim chirp loop

% =========================================================================
%% ── Step 4: CPI-Level Aggregation ───────────────────────────────────────
% =========================================================================

% Per-CPI interference ratio: fraction of victim chirps with any overlap.
if K == 0
    interference_probability = 0;
else
    interference_probability = sum(mask_chirp) / K;
end

end   % ARIFOMA_analytic_detector_batsinspired

% =========================================================================
%% ── Local Helper Functions ───────────────────────────────────────────────
% =========================================================================

% -------------------------------------------------------------------------
function segs = chirp_to_segments(t, f, mode)
% CHIRP_TO_SEGMENTS
%   Convert one chirp's (t, f) sample sequence into monotonic linear segments.
%
% PURPOSE
% -------
%   An FMCW chirp consists of one or more monotonic ramp segments separated
%   by direction changes (ramp → flyback → return).  This function detects
%   segment boundaries by sign changes in the instantaneous slope and fits
%   a linear model f(t) = a*t + b to each segment.
%
% MODES
%   'primary' — return only the main TX ramp segment (longest segment,
%               with a tie-break favouring the median slope direction).
%               Used for victim chirps and TX-only interferer modelling.
%   'all'     — return all monotonic segments (ramp + flyback + return).
%               Used when interferer flyback is included in the model.
%
% INPUTS
%   t     (double array) — sample times [s]
%   f     (double array) — instantaneous frequencies [Hz]
%   mode  (char)         — 'primary' or 'all'
%
% OUTPUT
%   segs  (struct array) — array of segment structs with fields:
%                            t0  (double) — segment start time [s]
%                            t1  (double) — segment end time [s]
%                            a   (double) — slope [Hz/s]
%                            b   (double) — intercept [Hz]
%                          Empty if fewer than 2 valid samples are available.

    t = t(:);  f = f(:);
    segs = struct('t0',{}, 't1',{}, 'a',{}, 'b',{});

    % ── Remove non-finite samples ─────────────────────────────────────────
    ok = isfinite(t) & isfinite(f);
    t  = t(ok);  f = f(ok);
    if numel(t) < 2,  return;  end

    % ── Sort by time and enforce unique timestamps ────────────────────────
    [t, ord] = sort(t, 'ascend');  f = f(ord);
    [tU, ia] = unique(t, 'stable');  fU = f(ia);
    t = tU;  f = fU;
    if numel(t) < 2,  return;  end

    % ── Compute instantaneous slopes ─────────────────────────────────────
    dt    = diff(t);
    df    = diff(f);
    valid = isfinite(dt) & isfinite(df) & (dt > 0);
    if nnz(valid) < 1,  return;  end

    slope = df(valid) ./ dt(valid);

    % Map valid slopes back to the full diff-index array.
    s        = zeros(size(df));
    s(valid) = slope;

    % ── Slope Sign Segmentation ───────────────────────────────────────────
    %   Small slopes below the numerical tolerance are treated as zero to
    %   prevent spurious segment boundaries from floating-point noise.

    tol_slope = 1e-12 * max(1, max(abs(slope)));
    sgn       = sign(s);
    sgn(abs(s) < tol_slope) = 0;

    % Propagate non-zero sign values forward then backward to fill zeros.
    sgnFilled = sgn;
    for i = 2:numel(sgnFilled)
        if sgnFilled(i) == 0,  sgnFilled(i) = sgnFilled(i-1);  end
    end
    for i = numel(sgnFilled)-1 : -1 : 1
        if sgnFilled(i) == 0,  sgnFilled(i) = sgnFilled(i+1);  end
    end
    if all(sgnFilled == 0),  sgnFilled(:) = 1;  end

    % ── Detect Segment Boundaries ─────────────────────────────────────────
    bDiff  = find(sgnFilled(1:end-1) ~= sgnFilled(2:end));
    bSamp  = bDiff + 1;   % sample index where each new segment begins

    starts = [1;      bSamp(:)];
    ends   = [bSamp(:); numel(t)];

    % ── Fit Linear Model to Each Segment ─────────────────────────────────
    tmp = struct('t0',{}, 't1',{}, 'a',{}, 'b',{}, 'dur',{}, 'dir',{});

    for k = 1:numel(starts)
        i0 = starts(k);
        i1 = ends(k);

        if (i1 - i0) < 1,  continue;  end   % need at least 2 points

        t0_seg = t(i0);  t1_seg = t(i1);
        f0_seg = f(i0);  f1_seg = f(i1);

        if ~(isfinite(t0_seg) && isfinite(t1_seg) && ...
             isfinite(f0_seg) && isfinite(f1_seg)) || (t1_seg <= t0_seg)
            continue;
        end

        a_seg = (f1_seg - f0_seg) / (t1_seg - t0_seg);
        b_seg = f0_seg - a_seg * t0_seg;

        tmp(end+1).t0  = t0_seg;          %#ok<AGROW>
        tmp(end).t1    = t1_seg;
        tmp(end).a     = a_seg;
        tmp(end).b     = b_seg;
        tmp(end).dur   = t1_seg - t0_seg;
        tmp(end).dir   = sign(a_seg);
    end

    if isempty(tmp),  return;  end

    % ── Mode: Return All Segments ─────────────────────────────────────────
    if strcmpi(mode, 'all')
        segs = rmfield(tmp, {'dur','dir'});
        return;
    end

    % ── Mode: Return Primary (Main TX Ramp) Segment ───────────────────────
    %   Selection heuristic:
    %     1. Score by segment duration (TX ramp is typically the longest).
    %     2. Tie-break by matching the median slope direction of the chirp.
    %   This correctly identifies the main ramp for positive-slope,
    %   negative-slope, and stepped FMCW waveforms.

    dirs   = [tmp.dir];
    durs   = [tmp.dur];

    dirMed = sign(median(slope(isfinite(slope))));
    if dirMed == 0,  dirMed = 1;  end

    dirMatch = (dirs == dirMed);
    score    = durs + 1e-6 * dirMatch;   % tiny direction-match bias
    [~, imax] = max(score);

    segs = rmfield(tmp(imax), {'dur','dir'});
end

% -------------------------------------------------------------------------
function J = abs_linear_interval(p, q, thr, tL, tR)
% ABS_LINEAR_INTERVAL
%   Solve |p*t + q| <= thr over the interval [tL, tR].
%
% PURPOSE
% -------
%   Given the linear IF frequency difference D(t) = p*t + q between a
%   victim and interferer chirp segment, returns the sub-interval of
%   [tL, tR] during which |D(t)| <= thr (i.e. the interferer falls within
%   the victim's LPF bandwidth).
%
%   The constraint |p*t + q| <= thr defines a convex set (intersection of
%   two half-spaces), so the result is always a single interval or empty.
%
% INPUTS
%   p, q  (double) — linear model coefficients: D(t) = p*t + q [Hz, Hz/s]
%   thr   (double) — bandwidth threshold (BADC) [Hz]
%   tL    (double) — interval left bound [s]
%   tR    (double) — interval right bound [s]
%
% OUTPUT
%   J  (double, 1×2 or 0×2) — [ts, te] if a valid overlap exists; else []

    p   = double(p);    q   = double(q);
    thr = double(thr);  tL  = double(tL);   tR = double(tR);

    J = zeros(0, 2);
    if tR <= tL,  return;  end

    if abs(p) < eps
        % ── Constant difference: D(t) = q ────────────────────────────────
        %   Overlap exists over the entire window if |q| <= thr.
        if abs(q) <= thr
            J = [tL, tR];
        end
        return;
    end

    % ── Linear difference: solve -thr <= p*t + q <= thr ──────────────────
    %   Roots of p*t + q = ±thr:
    t1  = (-thr - q) / p;
    t2  = ( thr - q) / p;
    low  = min(t1, t2);
    high = max(t1, t2);

    % Intersect solution set with [tL, tR].
    ts = max(tL, low);
    te = min(tR, high);

    if te > ts
        J = [ts, te];
    end
end

% -------------------------------------------------------------------------
function I = merge_intervals(Iin, tol)
% MERGE_INTERVALS
%   Merge overlapping or adjacent intervals into a minimal covering set.
%
% PURPOSE
% -------
%   When multiple interferer segments produce overlapping detection windows
%   for the same victim chirp, merging prevents double-counting in the
%   duration and centroid calculations.
%
% INPUTS
%   Iin  (double, N×2) — interval matrix with columns [start, end];
%                        rows need not be sorted or non-overlapping
%   tol  (double)      — adjacency tolerance [s]; intervals whose gap is
%                        less than or equal to tol are merged
%
% OUTPUT
%   I    (double, M×2) — merged interval matrix (M <= N), sorted by start time

    if isempty(Iin)
        I = Iin;
        return;
    end

    tol  = double(tol);
    Iin  = sortrows(Iin, 1);
    I    = Iin(1, :);

    for k = 2:size(Iin, 1)
        s = Iin(k, 1);
        e = Iin(k, 2);

        if s <= I(end, 2) + tol
            % Overlapping or adjacent: extend the current interval.
            I(end, 2) = max(I(end, 2), e);
        else
            % Disjoint: start a new interval.
            I(end+1, :) = [s, e];   %#ok<AGROW>
        end
    end
end