function [cornerData, pathByCorner, links_sym, radarCtx] = ARIFOMA_prepare_scene_views( ...
    out, vLoc_table, Param, corners, radarCtx, tolLinks, cfg)
%PREPARE_SCENE_VIEWS
% Build all tables/structures needed by do_tile:
%   - links + symmetric links
%   - per-corner victim + interferers (old format)
%   - per-sensor "perspective" and per-sensor "pairs" (new format)
%   - pathByCorner (used for distance extraction)
%   - sensorSeedLUT stored into radarCtx once
%
% Inputs
%   out         : frame struct from frame_to_mc() (must contain out.vehId)
%   vLoc_table  : scenario vLoc table
%   Param       : radar parameter table (VehID, Corner, SensorID, HasSensor, ...)
%   corners     : string array ["FL","FR","RL","RR"]
%   radarCtx    : your persistent radarCtx struct (will receive LUTs)
%   tolLinks    : time tolerance for link extraction (e.g. 0.002)
%
% Outputs
%   cornerData  : struct array, one element per corner
%   pathByCorner: your path range structure
%   links_sym   : symmetric link table (both directions)
%   radarCtx    : updated with sensorSeedLUT + waveformBySensor cache
 
% -------------------------------------------------------------
% A) Build links table for the current frame + per-corner path ranges
% -------------------------------------------------------------
if nargin < 7,  cfg = struct();  end
links       = ARIFOMA_build_interference_link_table_from_frame(out, vLoc_table, Param, tolLinks, cfg);
pathByCorner = ARIFOMA_build_interference_path_ranges_per_corner(links, out.vehId, corners);
 
% -------------------------------------------------------------
% B) Sensor lookup tables (stable scene-wide mapping)
% -------------------------------------------------------------
sensorLUT = table(Param.SensorID, Param.VehID, Param.Corner, ...
    'VariableNames', {'SensorID','VehID','Corner'});
% keep only links where ego vehicle is the victim
links     = enrich_links_with_ids(links, sensorLUT);
links_sym = make_links_symmetric(links);
 
% -------------------------------------------------------------
% C) Build sensorSeedLUT ONCE (scene-wide, not per-corner)
% -------------------------------------------------------------
radarCtx.sensorSeedLUT     = build_sensor_seed_lut(Param);
radarCtx.waveformBySensor  = struct();  % cache: waveform simulation results by SensorID
 
% -------------------------------------------------------------
% D) Build per-corner data (victim + interferers + per-sensor perspectives)
% -------------------------------------------------------------
cornerData = repmat(struct( ...
    'corner',        "", ...
    'victim',        [], ...
    'interferers',   struct([]), ...
    'perspective',   struct(), ...
    'pairsBySensor', struct() ), numel(corners), 1);
 
    for c = 1:numel(corners)
        corner = corners(c);
    
        % Victim for this corner on this vehicle
        Pv = Param(Param.VehID==out.vehId & Param.Corner==corner & Param.HasSensor, :);
        if isempty(Pv)
            cornerData(c).corner          = corner;
            cornerData(c).victim          = [];
            cornerData(c).interferers     = [];
            cornerData(c).interferersAll  = [];
            cornerData(c).perspective     = struct();
            cornerData(c).pairsBySensor   = struct();
            cornerData(c).pairs           = links_sym;   % stash symmetric links anyway
            continue;
        end
    
        % -----------------------------
        % 1) Victim-centric paths (this is what you already had)
        % -----------------------------
        Lc = links( double(links.VictimVehID)==double(out.vehId) & ...
            string(links.VictimCorner)==string(corner), :);
        sidList = unique(Lc.InterfererSensorID);
        sidList = sidList(isfinite(sidList));
        if isempty(sidList)
            PI = Param([],:);
        else
            PI = Param(ismember(Param.SensorID, sidList) & Param.HasSensor, :);
        end
    
        % Optional path stats per interferer (victim-centric)
        agg_victim = table();
        if ~isempty(Lc)
            assert(ismember('InterfererSensorID', Lc.Properties.VariableNames), ...
                   'links missing InterfererSensorID');
            if ~ismember('d1_m', Lc.Properties.VariableNames),  Lc.d1_m  = NaN(height(Lc),1); end
            if ~ismember('d2_m', Lc.Properties.VariableNames),  Lc.d2_m  = NaN(height(Lc),1); end
    
            [uSID, ~, G] = unique(Lc.InterfererSensorID);
            CountPaths = accumarray(G, 1, [numel(uSID), 1], @sum, 0);
            min_d1_m   = accumarray(G, Lc.d1_m, [numel(uSID), 1], @(x) min(x, [], 'omitnan'), NaN);
            min_d2_m   = accumarray(G, Lc.d2_m, [numel(uSID), 1], @(x) min(x, [], 'omitnan'), NaN);
            agg_victim = table(uSID, CountPaths, min_d1_m, min_d2_m, ...
                'VariableNames', {'InterfererSensorID','CountPaths','min_d1_m','min_d2_m'});
        end
    
        % Victim + primary interferers as in your original code
        victim = row2struct(Pv(1,:));
        if isempty(PI)
            interferers = struct([]);
        else
            interferers = arrayfun(@(i) param_row_to_seed(PI(i,:)), 1:height(PI));
        end
    
        if ~isempty(interferers) && ~isempty(agg_victim)
            for k = 1:numel(interferers)
                sid = interferers(k).SensorID;
                arow = agg_victim(agg_victim.InterfererSensorID==sid, :);
                if ~isempty(arow)
                    interferers(k).CountPaths = arow.CountPaths;
                    interferers(k).min_d1_m   = arow.min_d1_m;
                    interferers(k).min_d2_m   = arow.min_d2_m;
                else
                    interferers(k).CountPaths = 0;
                    interferers(k).min_d1_m   = NaN;
                    interferers(k).min_d2_m   = NaN;
                end
            end
        end
    
        % Store victim + "primary" interferers for compatibility
        cornerData(c).corner      = corner;
        cornerData(c).victim      = victim;
        cornerData(c).interferers = interferers;
    
        % =====================================================================
        % 2) Build true per-sensor perspectives from PATHS (links_sym)
        %    and a superset of interferers for waveform generation
        % =====================================================================
        
    pairsBySensor = struct();
    perspective   = struct();
 
    % --- 1) Victim perspective ---
    if ~isempty(victim) && isfield(victim,'SensorID') && isfinite(victim.SensorID)
        sid_V = double(victim.SensorID);
        keyV  = sprintf('S%d', sid_V);
 
        P_V = build_perspective_for_sid(sid_V, links);
 
        perspective.(keyV)   = P_V;
        % If you want full “pairs”, use your existing builder; otherwise just store agg:
        pairsBySensor.(keyV) = build_pairs_for_actor(sid_V, P_V.interfererSIDs, P_V.agg);
    end
 
    % --- 2) Each interferer’s own perspective ---
    for k = 1:numel(interferers)
        sid_I = double(interferers(k).SensorID);
        if ~isfinite(sid_I), continue; end
 
        keyI = sprintf('S%d', sid_I);
 
        P_I = build_perspective_for_sid(sid_I, links);
 
        perspective.(keyI)   = P_I;
        pairsBySensor.(keyI) = build_pairs_for_actor(sid_I, P_I.interfererSIDs, P_I.agg);
    end
 
        cornerData(c).pairsBySensor = pairsBySensor;
        cornerData(c).perspective   = perspective;
    end
end
 
function [links, unmatchedVictimKeys] = enrich_links_with_ids(links, sensorLUT)
%ENRICH_LINKS_WITH_IDS (robust)
% - VictimSensorID: lookup by (VictimVehID, VictimCorner) using normalized keys
% - InterfererVehID/Corner: lookup by InterfererSensorID
%
% Returns unmatchedVictimKeys for debugging.
 
% ---------- normalize types ----------
vVeh  = double(links.VictimVehID);
vCor  = upper(strtrim(string(links.VictimCorner)));
 
lutVeh = double(sensorLUT.VehID);
lutCor = upper(strtrim(string(sensorLUT.Corner)));
 
% ---------- build lookup keys (VehID|CORNER) ----------
keyLinks = compose("%d|%s", vVeh, vCor);
keyLUT   = compose("%d|%s", lutVeh, lutCor);
 
% ---------- VictimSensorID ----------
[tfV, locV] = ismember(keyLinks, keyLUT);
 
links.VictimSensorID = nan(height(links),1);
links.VictimSensorID(tfV) = double(sensorLUT.SensorID(locV(tfV)));
 
unmatchedVictimKeys = unique(keyLinks(~tfV));
 
% ---------- InterfererVehID/Corner from InterfererSensorID ----------
[tfI, locI] = ismember(double(links.InterfererSensorID), double(sensorLUT.SensorID));
 
links.InterfererVehID  = nan(height(links),1);
links.InterfererCorner = strings(height(links),1);
 
links.InterfererVehID(tfI)  = double(sensorLUT.VehID(locI(tfI)));
links.InterfererCorner(tfI) = string(sensorLUT.Corner(locI(tfI)));
end
 
 
function links_sym = make_links_symmetric(links)
%MAKE_LINKS_SYMMETRIC
% Duplicate link rows with victim/interferer swapped.
 
Lrev = links;
 
% swap victim fields
Lrev.VictimVehID    = links.InterfererVehID;
Lrev.VictimCorner   = links.InterfererCorner;
Lrev.VictimSensorID = links.InterfererSensorID;
 
% swap interferer fields
Lrev.InterfererVehID    = links.VictimVehID;
Lrev.InterfererCorner   = links.VictimCorner;
Lrev.InterfererSensorID = links.VictimSensorID;
 
% Path metrics should exist and be symmetric
need = {'d1_m','d2_m','d_LOS','d_REFL_total'};
for k = 1:numel(need)
    f = need{k};
    if ~ismember(f, links.Properties.VariableNames), links.(f) = nan(height(links),1); end
    if ~ismember(f, Lrev.Properties.VariableNames),  Lrev.(f)  = nan(height(Lrev),1);  end
end
 
links_sym = [links; Lrev];
end
function agg = aggregate_paths_by_interferer(Lc)
%AGGREGATE_PATHS_BY_INTERFERER
% Returns table: InterfererSensorID, CountPaths, min_d1_m, min_d2_m
 
if isempty(Lc)
    agg = table();
    return;
end
 
if ~ismember('d1_m', Lc.Properties.VariableNames), Lc.d1_m = nan(height(Lc),1); end
if ~ismember('d2_m', Lc.Properties.VariableNames), Lc.d2_m = nan(height(Lc),1); end
 
[uSID, ~, G] = unique(Lc.InterfererSensorID);
 
CountPaths = accumarray(G, 1, [numel(uSID),1], @sum, 0);
min_d1_m   = accumarray(G, Lc.d1_m, [numel(uSID),1], @(x) min(x,[],'omitnan'), NaN);
min_d2_m   = accumarray(G, Lc.d2_m, [numel(uSID),1], @(x) min(x,[],'omitnan'), NaN);
 
agg = table(uSID, CountPaths, min_d1_m, min_d2_m, ...
    'VariableNames', {'InterfererSensorID','CountPaths','min_d1_m','min_d2_m'});
end
function interferers = attach_agg_to_interferers(interferers, agg)
%ATTACH_AGG_TO_INTERFERERS
% Adds CountPaths, min_d1_m, min_d2_m to each interferer struct.
 
for k = 1:numel(interferers)
    sid = double(interferers(k).SensorID);
    arow = agg(agg.InterfererSensorID==sid, :);
 
    if ~isempty(arow)
        interferers(k).CountPaths = arow.CountPaths;
        interferers(k).min_d1_m   = arow.min_d1_m;
        interferers(k).min_d2_m   = arow.min_d2_m;
    else
        interferers(k).CountPaths = 0;
        interferers(k).min_d1_m   = NaN;
        interferers(k).min_d2_m   = NaN;
    end
end
end
function [perspective, pairsBySensor] = build_all_perspectives_for_corner(victim, interferers, links)
%BUILD_ALL_PERSPECTIVES_FOR_CORNER
% For each actor sensor (victim + interferers), build:
%   perspective.S{id}   = build_perspective_for_sid(...)
%   pairsBySensor.S{id} = build_pairs_for_actor(...)
 
perspective   = struct();
pairsBySensor = struct();
 
% Victim as an actor
sidV = double(victim.SensorID);
keyV = sprintf('S%d', sidV);
P_V  = build_perspective_for_sid(sidV, links);
perspective.(keyV)   = P_V;
pairsBySensor.(keyV) = build_pairs_for_actor(sidV, P_V.interfererSIDs, P_V.agg);
 
% Each interferer as an actor
for k = 1:numel(interferers)
    sidI = double(interferers(k).SensorID);
    keyI = sprintf('S%d', sidI);
 
    P_I = build_perspective_for_sid(sidI, links);
 
    perspective.(keyI)   = P_I;
    pairsBySensor.(keyI) = build_pairs_for_actor(sidI, P_I.interfererSIDs, P_I.agg);
end
end
function sensorSeedLUT = build_sensor_seed_lut(Param)
%BUILD_SENSOR_SEED_LUT
% Struct where each field S{id} contains param_row_to_seed(row)
 
sensorSeedLUT = struct();
 
for i = 1:height(Param)
    if ismember('HasSensor', Param.Properties.VariableNames) && ~Param.HasSensor(i)
        continue;
    end
    sid = double(Param.SensorID(i));
    if ~isfinite(sid), continue; end
 
    key = sprintf('S%d', sid);
    if ~isfield(sensorSeedLUT, key)
        sensorSeedLUT.(key) = param_row_to_seed(Param(i,:));
    end
end
end
 
function S = param_row_to_seed(pr)
    if height(pr) ~= 1
        error('param_row_to_seed: expects a single-row table.');
    end
    S.SensorID   = val(pr,'SensorID', NaN);
    S.VehID      = val(pr,'VehID', NaN);
    S.Corner     = char(string(val(pr,'Corner',"")));
    S.Scheme     = char(string(val(pr,'Scheme',"")));
    S.Waveform   = char(string(val(pr,'Waveform',"")));
    S.WaveformClean = char(string(val(pr,'WaveformClean',"NORMAL")));
    S.Dir        = val(pr,'Dir', +1);
    S.F0         = val(pr,'F0', NaN);
    S.Fc         = val(pr,'Fc', NaN);
    S.B          = val(pr,'B',  NaN);
    S.T          = val(pr,'Tch',NaN);
    S.Nch        = val(pr,'Nch',NaN);
    S.FreqIncrHz = val(pr,'FreqIncrHz',0);
    S.ChirpsPerStep = val(pr,'ChirpsPerStep', NaN);
    S.bandMin    = val(pr,'bandMin', NaN);
    S.bandMax    = val(pr,'bandMax', NaN);
    S.T_flyback_default = val(pr,'T_flyback_default', NaN);
    S.T_flyback_stepped = val(pr,'T_flyback_stepped', NaN);
    S.T_returnfly       = val(pr,'T_returnfly', NaN);
    S.Ptx_dBm = val(pr,'Ptx_dBm', NaN);
    S.Gt_dBi  = val(pr,'Gt_dBi',  NaN);
    S.Gr_dBi  = val(pr,'Gr_dBi',  NaN);
    S.Nseg=double(pr.Nseg);
    S.PermPattern=pr.PermPattern;
end
 
function x = val(t, name, default)
    if ismember(name, t.Properties.VariableNames)
        v = t.(name);
        if iscell(v) && ~isempty(v), v = v{1}; end
        if istable(v), v = v{1,1}; end
        if isstring(v) || ischar(v)
            x = v(1);
        else
            x = v(1);
        end
    else
        x = default;
    end
end
function P = build_perspective_for_sid(sidVictim, links)
% build_perspective_for_sid  Return how a given sensor ID "sees" others,
% based purely on the links table.
%
% P.victimSID      : this sensor ID
% P.interfererSIDs : vector of SensorIDs that interfere with this sensor
% P.agg            : table with CountPaths, min_d1_m, min_d2_m per interferer
 
    % All rows where this sensor is victim (from any vehicle/corner)
    L = links(links.VictimSensorID == sidVictim & isfinite(links.InterfererSensorID), :);
 
    P = struct();
    P.victimSID      = sidVictim;
    P.interfererSIDs = [];
    P.agg            = table();
 
    if isempty(L)
        return;
    end
 
    [uSID, ~, G] = unique(L.InterfererSensorID);
 
    CountPaths = accumarray(G, 1, [numel(uSID), 1], @sum, 0);
 
    % Handle distance fields if present
    has_d1 = ismember('d1_m', L.Properties.VariableNames);
    has_d2 = ismember('d2_m', L.Properties.VariableNames);
 
    if has_d1 && has_d2
        min_d1_m = accumarray(G, L.d1_m, [numel(uSID), 1], ...
                              @(x) min(x, [], 'omitnan'), NaN);
        min_d2_m = accumarray(G, L.d2_m, [numel(uSID), 1], ...
                              @(x) min(x, [], 'omitnan'), NaN);
        agg = table(uSID, CountPaths, min_d1_m, min_d2_m, ...
            'VariableNames', {'InterfererSensorID','CountPaths','min_d1_m','min_d2_m'});
    else
        agg = table(uSID, CountPaths, ...
            'VariableNames', {'InterfererSensorID','CountPaths'});
    end
 
    P.interfererSIDs = uSID;
    P.agg            = agg;
end
function pairs = build_pairs_for_actor(actorSID, othersSID, agg)
% Build an ego-centric geometry map: "if ACTOR were the victim".
% 'agg' is your table with columns: InterfererSensorID, min_d1_m, min_d2_m.
    key = @(a,b) sprintf('S%d_S%d', round(double(a)), round(double(b)));
    pairs = struct();
    if isempty(othersSID), return; end
 
    for k = 1:numel(othersSID)
        sid = othersSID(k);
        d1 = NaN; d2 = NaN;
        if istable(agg) && ~isempty(agg)
            r = agg(agg.InterfererSensorID == sid, :);
            if ~isempty(r)
                if ismember('min_d1_m', r.Properties.VariableNames), d1 = r.min_d1_m(1); end
                if ismember('min_d2_m', r.Properties.VariableNames), d2 = r.min_d2_m(1); end
            end
        end
        rec = struct('A',actorSID,'B',sid,'d_LOS',double(d1),'d_REFL_total',double(d1+d2));
        pairs.(key(actorSID,sid)) = rec;   % directed: actor -> other
        pairs.(key(sid,actorSID)) = rec;   % keep symmetric for convenience
    end
end
 
function S = row2struct(pr)
    S = struct();
    S.VehID    = double(pr.VehID);
    S.Corner   = string(pr.Corner);
    S.SensorID = double(pr.SensorID);
    S.Scheme   = string(pr.Scheme);
    S.Waveform = string(pr.Waveform);
    S.WaveformClean = string(pr.WaveformClean);
    S.Dir      = double(pr.Dir);
    S.B        = double(pr.B);
    S.T      = double(pr.Tch);
    S.Nch      = double(pr.Nch);
    S.ChirpsPerStep = double(pr.ChirpsPerStep);
    S.FreqIncrHz    = double(pr.FreqIncrHz);
    S.F0       = double(pr.F0);
    S.Fc       = double(pr.Fc);
    S.bandMin  = double(pr.bandMin);
    S.bandMax  = double(pr.bandMax);
    S.T_flyback_default = double(pr.T_flyback_default);
    S.T_flyback_stepped = double(pr.T_flyback_stepped);
    S.T_returnfly       = double(pr.T_returnfly);
    S.Nseg=double(pr.Nseg);
    S.PermPattern=pr.PermPattern;
end