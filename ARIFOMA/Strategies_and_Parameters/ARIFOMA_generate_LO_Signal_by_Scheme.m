function [t_all, f_all, idx_all] = ARIFOMA_generate_LO_Signal_by_Scheme( ...
    scheme, waveform, f0, B, T, N, Tfly_default, Tfly_stepped, Treturn, ...
    chirpsPerStep, fIncr, dir, stepDir, bandMin, bandMax, Nseg, permPattern)

wfU = upper(strtrim(waveform));
isStepped = contains(wfU,'STEPPED');

% If you want: enforce dir from waveform string (optional safety)
if contains(wfU,'NEGATIVE-SLOPE')
    dir = -1;
elseif contains(wfU,'POSITIVE-SLOPE')
    dir = +1;
end

genStatic  = @(fstart) generateChirpSequenceNormal(fstart, B, T, N, Tfly_default, dir);

if isStepped
    genStepped = @() generateChirpSequenceStepped( ...
        f0, B, T, N, Tfly_stepped, Treturn, chirpsPerStep, fIncr, dir, stepDir, bandMin, bandMax);
end

switch upper(strtrim(scheme))

    case 'CONVENTIONAL'
        if isStepped
            [t_all,f_all,idx_all] = genStepped();
        else
            [t_all,f_all,idx_all] = genStatic(clampF0(f0,B,dir,bandMin,bandMax));
        end
    case {'BATSINSPIRED_FH_PER_CPI'}
        if isStepped
            [t_all,f_all,idx_all] = genStepped(clampF0(f0,B,dir,bandMin,bandMax));
        else
            [t_all,f_all,idx_all] = genStatic(clampF0(f0,B,dir,bandMin,bandMax));
        end
    case {'RANDOM_FH_PER_CPI'}
        if isStepped
            f0rand = randomF0ForStepped( ...
                B, dir, bandMin, bandMax, chirpsPerStep, fIncr, stepDir);
            [t_all,f_all,idx_all] = genStepped(f0rand);
        else
            f0rand = randomF0ForNormal(B, dir, bandMin, bandMax);
            [t_all,f_all,idx_all] = genStatic(f0rand);
        end

    case {'FAR'}
        if isStepped
            error('FAR not built for STEPPED');
        end
        [t_all,f_all,idx_all] = generateChirpSequenceRandPerChirp(B,T,N,Tfly_default,dir,bandMin,bandMax);

    case {'SEGMENTED'}
        sps = 32;
        [t_all, f_all, idx_all] = generateChirpSequenceSegmented( ...
            clampF0(f0, B, dir, bandMin, bandMax), B, T, N, Nseg, Tfly_default, ...
            dir, bandMin, bandMax, ...
            'permPattern', permPattern{1}, 'samplesPerSegment', sps);

        % sanitation
        ok = isfinite(t_all) & isfinite(f_all);
        t_all = t_all(ok); f_all = f_all(ok); idx_all = idx_all(ok);
        [t_all, ord] = sort(t_all, 'ascend');
        f_all = f_all(ord); idx_all = idx_all(ord);
        [tU, ia] = unique(t_all, 'stable');
        t_all = tU; f_all = f_all(ia); idx_all = idx_all(ia);

    otherwise
        error('Unknown scheme: %s', scheme);
end
end


function [t_all,f_all,idx_all] = generateChirpSequenceNormal(f0,B,T,N,T_fly,dir)
Nsamp = 100;
t_all=[]; f_all=[]; idx_all=[];
for k=1:N
    f_start=f0;
    f_end=f0 + dir*B;
    t_chirp=linspace(0,T,Nsamp);
    f_chirp=linspace(f_start,f_end,Nsamp);
    t_off=(k-1)*(T+T_fly);
    t_chirp=t_chirp+t_off;
    t_all=[t_all,t_chirp];
    f_all=[f_all,f_chirp];
    idx_all=[idx_all, k*ones(1,numel(t_chirp))];
end
end


function [t_all,f_all,idx_all] = generateChirpSequenceRandPerChirp(B,T,N,Tfly,dir,bandMin,bandMax)
Nsamp = 100;
t_all=[]; f_all=[]; idx_all=[];
for k=1:N
    if dir==1
        f0k = bandMin + (bandMax-bandMin-B)*rand;
        f_start=f0k; f_end=f0k+B;
    else
        f0k = bandMin + B + (bandMax-bandMin-B)*rand;
        f_start=f0k; f_end=f0k-B;
    end
    t_chirp=linspace(0,T,Nsamp);
    f_chirp=linspace(f_start,f_end,Nsamp);
    t_off=(k-1)*(T+Tfly);
    t_chirp=t_chirp+t_off;
    t_all=[t_all,t_chirp];
    f_all=[f_all,f_chirp];
    idx_all=[idx_all, k*ones(1,numel(t_chirp))];
end
end


function f0c = clampF0(f0,B,dir,bandMin,bandMax)
if dir==1
    f0c = min(max(f0, bandMin), bandMax - B);
else
    f0c = min(max(f0, bandMin + B), bandMax);
end
end
function [t_all,f_all,idx_all] = generateChirpSequenceStepped( ...
    f0,B,T,N,T_fly,T_ret,chirpsPerStep,fIncr,dir,stepDir,bandMin,bandMax)

Nsamp = 100;

t_all = [];
f_all = [];
idx_all = [];

t0 = 0;  % running time pointer

for k = 1:N
    % step index repeats 0..(chirpsPerStep-1), then resets
    stepIdx = mod(k-1, chirpsPerStep);

    % start frequency increments each chirp within the cycle
    f_step  = f0 + stepDir*stepIdx*fIncr;

    % clamp so the entire chirp stays inside the band
    f_start = clampF0(f_step, B, dir, bandMin, bandMax);
    f_end   = f_start + dir*B;

    % samples
    t_chirp = linspace(0, T, Nsamp) + t0;
    f_chirp = linspace(f_start, f_end, Nsamp);

    t_all   = [t_all,  t_chirp];
    f_all   = [f_all,  f_chirp];
    idx_all = [idx_all, k*ones(1,Nsamp)];

    % advance time (continuous)
    t0 = t0 + T + T_fly;

    % after finishing a full cycle of chirpsPerStep, add return time once
    if mod(k, chirpsPerStep) == 0
        t0 = t0 + T_ret;
    end
end
end
function f0r = randomF0ForStepped(B, dir, bandMin, bandMax, chirpsPerStep, fIncr, stepDir)

nSteps = max(1, round(chirpsPerStep));
df = abs(fIncr);

% valid start-frequency interval for one chirp
if dir >= 0
    fmin = bandMin;
    fmax = bandMax - B;
else
    fmin = bandMin + B;
    fmax = bandMax;
end

% widen constraint depending on step direction so all steps fit
if stepDir >= 0
    lo = fmin;
    hi = fmax - (nSteps - 1)*df;
else
    lo = fmin + (nSteps - 1)*df;
    hi = fmax;
end

if hi < lo
    error('No valid random f0 range for stepped waveform. Reduce B, fIncr, or chirpsPerStep.');
end

if hi == lo
    f0r = lo;
else
    f0r = lo + rand*(hi - lo);
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
    error('No valid random f0 range: B exceeds available band.');
end

if hi == lo
    f0r = lo;
else
    f0r = lo + rand*(hi - lo);
end
end

function [t_all, f_all, idx_all] = generateChirpSequenceSegmented( ...
    f0, B, T, N, Nseg, T_fly, dir, bandMin, bandMax, varargin)

% Generates a segmented FMCW chirp sequence.
%
% Each chirp:
%   - spans total bandwidth B
%   - lasts total duration T
%   - is divided into Nseg small chirps
%   - those small chirps are transmitted in permuted order
%
% Inputs:
%   f0, B, T, N, Nseg, T_fly, dir, bandMin, bandMax
%
% Optional name-value pairs:
%   'permPattern'       : permutation of 1:Nseg
%   'samplesPerSegment' : samples used per small chirp (default 32)
%
% Outputs:
%   t_all, f_all, idx_all

    p = inputParser;
    addParameter(p, 'permPattern', 1:Nseg);
    addParameter(p, 'samplesPerSegment', 32);
    parse(p, varargin{:});

    permPattern = p.Results.permPattern;
    sps         = p.Results.samplesPerSegment;

    % make sure permutation is a row vector
    permPattern = permPattern(:).';

    % fallback if invalid
    if numel(permPattern) ~= Nseg || ~isequal(sort(permPattern), 1:Nseg)
        permPattern = 1:Nseg;
    end

    % clamp global chirp start frequency
    f0 = clampF0(f0, B, dir, bandMin, bandMax);

    % per-segment bandwidth and duration
    Bseg = B / Nseg;
    Tseg = T / Nseg;

    t_all   = [];
    f_all   = [];
    idx_all = [];

    for k = 1:N

        t_chirp0 = (k-1) * (T + T_fly);

        for q = 1:Nseg

            % which original segment is transmitted at this time slot?
            seg = permPattern(q);

            % segment time location inside the current chirp
            t_seg0 = t_chirp0 + (q-1) * Tseg;

            % sample times in this segment
            % use endpoint-excluded sampling to avoid duplicate times
            tau = (0:sps-1) / sps * Tseg;
            t_seg = t_seg0 + tau;

            % original segment frequency range before permutation
            if dir >= 0
                % positive-slope chirp: f0 -> f0 + B
                f_seg_start = f0 + (seg-1)*Bseg;
                f_seg_end   = f_seg_start + Bseg;
            else
                % negative-slope chirp: f0 -> f0 - B
                f_seg_start = f0 - (seg-1)*Bseg;
                f_seg_end   = f_seg_start - Bseg;
            end

            % local mini-chirp
            f_seg = linspace(f_seg_start, f_seg_end, sps);

            t_all   = [t_all, t_seg];
            f_all   = [f_all, f_seg];
            idx_all = [idx_all, k * ones(1, sps)];
        end
    end
end

