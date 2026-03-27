function frame = ARIFOMA_frame_extract_mc(vLoc, intf, Param, t_query, opt)
% ARIFOMA_frame_extract_mc
% Extract a single time snapshot for victim-centric MC processing (no plotting).
%
% Key behavior:
%   - Role filter: victim/source/either
%   - Optional X-window filter (WinX)
%   - Sensor-existence filter using Param.HasSensor
%   - Path length filter + max-path downselect
%
% New optional mode:
%   opt.SingleSensorMode (true/false):
%       If true, all vehicles are assumed to have exactly one sensor corner:
%       opt.SingleCornerLabel (default "FR"). Corner inference from (x,y) is skipped.
%
% Required opt fields (as before):
%   VehicleId, Role, Corner, Tol, MaxPaths, WinX, ShowAllVehicles, MaxPathLen
%
% Output:
%   frame struct with R_draw (filtered paths) + victimRows + interfererRows.

% ---------------- basic constants ----------------
allCorners = ["FL","FR","RR","RL"];  % only used when not in SingleSensorMode

if ~isfield(opt,'SingleSensorMode'),  opt.SingleSensorMode = false; end
if ~isfield(opt,'SingleCornerLabel'), opt.SingleCornerLabel = "FR"; end

% ---------------- 1) snap time and slice ----------------
[V, M, R, ~, ~, t0] = slice_wilab_at_time(vLoc, intf, t_query, opt.Tol);
vRows = normalize_vehicle_rows(V, M, t0);

% ---------------- 2) choose victim vehicle ----------------
vehId = opt.VehicleId;
if isempty(vehId)
    vehId = auto_pick_vehicle_in_gate(vRows, [0 100], t0);
end
vrVict = vrow_by_vid(vRows, vehId);
if isempty(vrVict)
    error('Vehicle %g not present at t≈%.3f.', vehId, t0);
end

% ---------------- 3) role filter ----------------
role  = lower(string(opt.Role));
Rrole = filter_paths_by_role(R, vehId, role);

% ---------------- 4) optional X-window around victim ----------------
XL = [];
if ~isempty(opt.WinX)
    cx = double(vrVict.V_x(1));
    XL = [cx - opt.WinX/2, cx + opt.WinX/2];
end

% ---------------- 5) victim sensor params (HasSensor only) ----------------
Pv_all = get_victim_param_rows(Param, vehId);

% If user wants single-sensor mode, force victim corner to that one label.
if opt.SingleSensorMode
    wantCorners = string(opt.SingleCornerLabel);
else
    wantCorners = choose_victim_corners(opt.Corner, allCorners, vrVict, Rrole, vehId, Pv_all);
end

victimRows = build_victim_rows(wantCorners, Pv_all, vehId, allCorners, opt.SingleSensorMode);

% ---------------- 6) start R_draw from role-filtered ----------------
R_draw = Rrole;

% ---------------- 7) apply window filter (vehicles + paths) ----------------
[Vin, R_draw] = apply_window_filter(vRows, R_draw, XL);

if ~opt.ShowAllVehicles
    Vin = Vin(Vin.Vid == vehId, :);
end

% ---------------- 8) RX corner filtering (only if NOT single sensor mode) ----------------
% In SingleSensorMode, every interferedSensor is treated as the single corner,
% so we do not filter by geometry.
if ~opt.SingleSensorMode
    R_draw = filter_rx_paths_by_corner(R_draw, vehId, vrVict, wantCorners, allCorners, opt.Corner);
end

% ---------------- 9) TX sensor existence filtering (Param.HasSensor) ----------------
% In SingleSensorMode we simply require that each interfering vehicle has the single corner sensor.
R_draw = filter_tx_paths_by_sensor_existence(R_draw, vRows, Param, allCorners, opt.SingleSensorMode, opt.SingleCornerLabel);

% ---------------- 10) Path length filter ----------------
R_draw = filter_paths_by_length(R_draw, opt.MaxPathLen);

% ---------------- 11) Downselect if too many paths (keep SHORTEST first) ----------------
R_draw = downselect_paths_if_needed(R_draw, opt.MaxPaths);

% ---------------- 12) Extract interferer param rows + seeds ----------------
[Irows, Iseeds] = build_interferer_rows_and_seeds(R_draw, vRows, Param, allCorners, opt.SingleSensorMode, opt.SingleCornerLabel);

% ---------------- 13) Pack output ----------------
frame = struct();
frame.t0          = t0;
frame.vehId       = vehId;
frame.role        = role;
frame.wantCorners = wantCorners;

frame.vRows  = vRows;
frame.vrVict = vrVict;
frame.Vin    = Vin;
frame.XL     = XL;

frame.R_all  = R;
frame.R_role = Rrole;
frame.R_draw = R_draw;

frame.victimRows     = victimRows;
frame.interfererRows = Irows;

frame.victimSeeds     = rowfun_to_seeds(victimRows);
frame.interfererSeeds = Iseeds;

end

%% ============================= Helpers =============================

function vRows = normalize_vehicle_rows(V, M, t0)
% Normalize vLoc slice into a table with WiLab corner schema.
if istable(V) && ~isempty(V)
    vRows = V;
    return;
end
if ~isempty(M)
    vRows = array2table(M(:,1:13), ...
        'VariableNames', {'time','Vid','V_x','V_y','V_speed', ...
        'corner_p1_x','corner_p1_y','corner_p2_x','corner_p2_y', ...
        'corner_p3_x','corner_p3_y','corner_p4_x','corner_p4_y'});
    return;
end
error('No vehicles at t≈%.3f.', t0);
end

function vehId = auto_pick_vehicle_in_gate(vRows, gateX, t0)
% Pick a victim inside a spatial gate (fallback).
xmin = gateX(1); xmax = gateX(2);
cand = vRows(vRows.V_x>=xmin & vRows.V_x<=xmax, :);
if isempty(cand)
    error('No vehicle in [%.1f, %.1f] m at t≈%.3f. Provide opt.VehicleId.', xmin, xmax, t0);
end
[~,ix] = min(cand.V_x);
vehId = double(cand.Vid(ix));
end

function Rrole = filter_paths_by_role(R, vehId, role)
% Keep only interference rows relevant to the chosen role.
switch role
    case "victim"
        Rrole = R(R.interferedId==vehId, :);
    case "source"
        Rrole = R(R.interferingId==vehId, :);
    otherwise
        Rrole = R(R.interferedId==vehId | R.interferingId==vehId, :);
end
end

function Pv_all = get_victim_param_rows(Param, vehId)
% Always prefer sensor-enabled rows.
if ismember('HasSensor', Param.Properties.VariableNames)
    Pv_all = Param(Param.VehID==vehId & Param.HasSensor, :);
else
    Pv_all = Param(Param.VehID==vehId, :);
end
end

function wantCorners = choose_victim_corners(cornerIn, allCorners, vrVict, Rrole, vehId, Pv_all)
% Corner selection without "front-only" concept.
% Supports: "all" | "auto" | explicit label | numeric 1..4

if isnumeric(cornerIn) && isscalar(cornerIn)
    cStr = allCorners(double(max(1,min(4,cornerIn))));
else
    cStr = upper(string(cornerIn));
end

switch cStr
    case "AUTO"
        % pick busiest victim sensorId if available
        if ismember('interferedSensorId', Rrole.Properties.VariableNames)
            ids = Rrole.interferedSensorId(Rrole.interferedId==vehId);
            ids = double(ids(isfinite(ids)));
        else
            ids = [];
        end
        if ~isempty(ids)
            u = unique(ids);
            cnt = arrayfun(@(x) sum(ids==x), u);
            sid = u(argmax(cnt));
            r1 = Rrole(find(double(Rrole.interferedSensorId)==sid,1,'first'),:);
            kc = corner_index_from_xy(vrVict, r1.interferedSensor_x, r1.interferedSensor_y);
            wantCorners = allCorners(kc);
        else
            % fallback: use Param corners if available, else FL
            if ~isempty(Pv_all) && ismember('Corner', Pv_all.Properties.VariableNames)
                wantCorners = string(Pv_all.Corner(1));
            else
                wantCorners = allCorners(1);
            end
        end

    case "ALL"
        if ~isempty(Pv_all) && ismember('Corner', Pv_all.Properties.VariableNames)
            wantCorners = unique(string(Pv_all.Corner),'stable');
        else
            wantCorners = allCorners;
        end

    otherwise
        wantCorners = string(cStr);
end
end

function victimRows = build_victim_rows(wantCorners, Pv_all, vehId, allCorners, singleMode)
% Build victim parameter rows for the selected corners.
%
% In singleMode, we only keep 1 corner; if Param doesn't contain it, we error
% (because single-sensor mode should be consistent with Param generation).

victimRows = table();
for cc = 1:numel(wantCorners)
    thisC = wantCorners(cc);

    Pv = Pv_all;
    if ~isempty(Pv) && ismember('Corner', Pv.Properties.VariableNames)
        Pv = Pv(Pv.Corner==thisC, :);
    end

    if isempty(Pv)
        if singleMode
            error('Victim %g has no Param row for single-sensor corner "%s".', vehId, thisC);
        else
            % minimal fallback only for development
            Pv = synth_param_row(vehId, index_from_label(thisC, allCorners), "LRR");
            Pv.Corner = thisC;
        end
    end

    victimRows = [victimRows; Pv(1,:)]; %#ok<AGROW>
end
end

function R_draw = filter_rx_paths_by_corner(R_draw, vehId, vrVict, wantCorners, allCorners, cornerOpt)
% Filter victim-side (RX) paths by selected corners (geometry-based).
needFilter = any(R_draw.interferedId==vehId) && ~strcmpi(string(cornerOpt),'all');
if ~needFilter, return; end

keepRX = false(height(R_draw),1);
for i = 1:height(R_draw)
    if R_draw.interferedId(i) ~= vehId
        continue;
    end
    k   = corner_index_from_xy(vrVict, R_draw.interferedSensor_x(i), R_draw.interferedSensor_y(i));
    lab = allCorners(k);
    keepRX(i) = any(lab == wantCorners);
end
R_draw = R_draw(keepRX,:);
end

function R_draw = filter_tx_paths_by_sensor_existence(R_draw, vRows, Param, allCorners, singleMode, singleCorner)
% Keep only interference rows where the interferer has an active sensor.
%
% In singleMode, every vehicle is assumed to have only 'singleCorner',
% so we only check that (VehID, singleCorner) exists and HasSensor==true.

if isempty(R_draw), return; end
if ~ismember('HasSensor', Param.Properties.VariableNames)
    return; % can't enforce
end

keepTX = true(height(R_draw),1);

for i = 1:height(R_draw)
    txVeh = double(R_draw.interferingId(i));
    if ~isfinite(txVeh), continue; end

    if singleMode
        lab = string(singleCorner);
    else
        vrTx  = vrow_by_vid(vRows, txVeh);
        if isempty(vrTx), continue; end
        ktx  = corner_index_from_xy(vrTx, R_draw.interferingSensor_x(i), R_draw.interferingSensor_y(i));
        lab  = allCorners(ktx);
    end

    Ptx = Param(Param.VehID==txVeh & Param.Corner==lab, :);
    if isempty(Ptx) || ~any(Ptx.HasSensor)
        keepTX(i) = false;
    end
end

R_draw = R_draw(keepTX,:);
end

function [Vin, R_draw] = apply_window_filter(vRows, R_draw, XL)
% Apply an X-window to vehicles and path endpoints (victim or interferer sensor x).
Vin = vRows;
if isempty(XL), return; end

Vin = Vin(Vin.V_x >= XL(1) & Vin.V_x <= XL(2), :);
if isempty(R_draw), return; end

keepI = (R_draw.interferedSensor_x >= XL(1) & R_draw.interferedSensor_x <= XL(2)) | ...
        (R_draw.interferingSensor_x>= XL(1) & R_draw.interferingSensor_x<= XL(2));
R_draw = R_draw(keepI,:);
end

function R_draw = filter_paths_by_length(R_draw, maxLen)
% Drop paths whose total path length exceeds maxLen.
if isempty(R_draw), return; end
if ~all(ismember({'lengthPathFirst','lengthPathSecond'}, R_draw.Properties.VariableNames))
    return;
end

L1 = double(R_draw.lengthPathFirst);
L2 = double(R_draw.lengthPathSecond);
L2(~isfinite(L2)) = 0;

totalLen = L1 + max(0, L2);
R_draw = R_draw(totalLen <= maxLen, :);
end

function R_draw = downselect_paths_if_needed(R_draw, maxPaths)
% Keep only the most relevant paths if too many exist.
% Relevance proxy here: shortest total path length first.
if height(R_draw) <= maxPaths, return; end

L1 = double(R_draw.lengthPathFirst);
L2 = double(R_draw.lengthPathSecond);
L2(~isfinite(L2)) = 0;
totalLen = L1 + max(0, L2);

[~,ord] = sort(totalLen,'ascend');   % shortest first (usually strongest)
R_draw = R_draw(ord(1:maxPaths), :);
end

function [Irows, Iseeds] = build_interferer_rows_and_seeds(R_draw, vRows, Param, allCorners, singleMode, singleCorner)
% Build unique interferer parameter rows from R_draw.
Irows  = table();
Iseeds = {};

if isempty(R_draw), return; end

for i = 1:height(R_draw)
    txVeh = double(R_draw.interferingId(i));
    if ~isfinite(txVeh), continue; end

    if singleMode
        labTx = string(singleCorner);
    else
        vrTx = vrow_by_vid(vRows, txVeh);
        if isempty(vrTx), continue; end
        ktx   = corner_index_from_xy(vrTx, R_draw.interferingSensor_x(i), R_draw.interferingSensor_y(i));
        labTx = allCorners(ktx);
    end

    Ptx = Param(Param.VehID==txVeh & Param.Corner==labTx, :);
    if isempty(Ptx)
        if singleMode
            % In single-sensor mode this should not happen if Param is consistent.
            continue;
        else
            Ptx = synth_param_row(txVeh, index_from_label(labTx, allCorners), "LRR");
            Ptx.Corner = labTx;
        end
    else
        if ismember('HasSensor', Ptx.Properties.VariableNames) && ~any(Ptx.HasSensor)
            continue;
        end
        Ptx = Ptx(1,:);
    end

    Irows = [Irows; Ptx]; %#ok<AGROW>
    Iseeds{end+1,1} = param_row_to_seed(Ptx); %#ok<AGROW>
end

% Deduplicate interferers by SensorID if available; else (VehID,Corner)
if ~isempty(Irows)
    if ismember('SensorID', Irows.Properties.VariableNames)
        [~,iu] = unique(Irows.SensorID,'stable');
    else
        [~,iu] = unique(string(Irows.VehID) + "_" + string(Irows.Corner), 'stable');
    end
    Irows  = Irows(iu,:);
    Iseeds = Iseeds(iu);
end
end

function seeds = rowfun_to_seeds(T)
seeds = cell(0,1);
for i=1:height(T)
    seeds{end+1,1} = param_row_to_seed(T(i,:)); %#ok<AGROW>
end
end

function S = param_row_to_seed(PR)
% Minimal seed struct (keep only what your PRI engine needs).
S = struct();
S.SensorID      = getIf(PR,'SensorID',NaN);
S.VehID         = PR.VehID;
S.Corner        = string(getIf(PR,'Corner',"FL"));
S.Scheme        = string(getIf(PR,'Scheme',"STATIC"));
S.Waveform      = string(getIf(PR,'Waveform',"NORMAL"));
S.WaveformClean = string(getIf(PR,'WaveformClean',"NORMAL"));
S.Dir           = getIf(PR,'Dir',+1);
S.f0            = getIf(PR,'F0',77e9);
S.B             = getIf(PR,'B',400e6);
S.T             = getIf(PR,'Tch',30e-6);
S.Nch           = getIf(PR,'Nch',16);
S.ChirpsPerStep = getIf(PR,'ChirpsPerStep',NaN);
S.FreqIncrHz    = getIf(PR,'FreqIncrHz',0);
S.bandMin       = getIf(PR,'bandMin',77e9);
S.bandMax       = getIf(PR,'bandMax',81e9);
S.T_fb_def      = getIf(PR,'T_flyback_default',3e-6);
S.T_fb_step     = getIf(PR,'T_flyback_stepped',0.5e-6);
S.T_return      = getIf(PR,'T_returnfly',4e-6);
end

function v = getIf(T, name, def)
if ismember(name, T.Properties.VariableNames)
    v = T.(name)(1);
else
    v = def;
end
end

function [V,M,R,rowsV,rowsI,t0] = slice_wilab_at_time(vLoc, intf, t_query, tol)
% Snap to nearest timestamp and slice both vLoc and intf.
if nargin<4 || isempty(tol), tol = 0.05; end

isVTable = istable(vLoc);
if isVTable
    tv = unique(double(vLoc.time));
else
    tv = unique(double(vLoc(:,1)));
end

ti   = unique(double(intf.time));
allT = unique([tv(:); ti(:)]);
assert(~isempty(allT),'No timestamps present.');

[~,k] = min(abs(allT - t_query));
t0 = allT(k);

if isVTable
    rowsV = abs(double(vLoc.time) - t0) <= tol;
    V = vLoc(rowsV,:);
    M = [];
else
    rowsV = abs(double(vLoc(:,1)) - t0) <= tol;
    M = vLoc(rowsV,:);
    V = table();
end

rowsI = abs(double(intf.time) - t0) <= tol;
R = intf(rowsI,:);
end

function vr = vrow_by_vid(vRows, vehId)
ix = find(double(vRows.Vid)==vehId, 1, 'first');
if isempty(ix), vr = []; else, vr = vRows(ix,:); end
end

function idx = argmax(v)
[~,idx] = max(v);
end

function k = corner_index_from_xy(vrow, sx, sy)
% Infer nearest corner index based on (sx,sy) compared to vrow corner points.
C = [vrow.corner_p1_x vrow.corner_p1_y; ...
     vrow.corner_p2_x vrow.corner_p2_y; ...
     vrow.corner_p3_x vrow.corner_p3_y; ...
     vrow.corner_p4_x vrow.corner_p4_y];
d = hypot(double(C(:,1))-double(sx), double(C(:,2))-double(sy));
[~,k] = min(d);
end

function n = index_from_label(lbl, allCorners)
lbl = upper(string(lbl));
m = find(allCorners==lbl,1);
if isempty(m), m = 1; end
n = m;
end

function T = synth_param_row(vehId, cornerIdx, scenario)
% Development-only fallback. For production, replace with error().
cornerList = ["FL","FR","RR","RL"];
S = struct('VehID',vehId,'Corner',cornerList(min(4,max(1,cornerIdx))), ...
    'HasSensor',true,'Scenario',string(scenario), ...
    'Scheme',"STATIC",'Waveform',"NORMAL",'WaveformClean',"NORMAL",'Dir',+1, ...
    'B',400e6,'Tch',30e-6,'Nch',16,'ChirpsPerStep',NaN,'FreqIncrHz',0, ...
    'F0',77e9,'Fc',77e9+200e6,'T_flyback_default',3e-6,'T_flyback_stepped',0.5e-6, ...
    'T_returnfly',4e-6,'bandMin',77e9,'bandMax',81e9);
T = struct2table(S);
end