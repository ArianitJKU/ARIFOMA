function extract_wilab_db_to_mat(baseDir)
% EXTRACT_WILAB_DB_TO_MAT
%   Convert WiLabVISim SQLite output databases to MATLAB .mat files.
%
% PURPOSE
% -------
%   WiLabVISim stores simulation results in SQLite databases
%   (output_database_1.db). This utility reads the two tables that ARIFOMA
%   requires — vLocation and firstReflectLog — converts them to MATLAB
%   tables and numeric matrices, and saves everything to a .mat file
%   (wilab_export.mat) in the same folder as the database.
%
%   The exported .mat file can be loaded directly by main_sim.m via
%   cfg.inputMat.
%
% USAGE
% -----
%   extract_wilab_db_to_mat('<path-to-run-folder>')
%
%   If baseDir is omitted, the function defaults to:
%     ./Output/antenna4/dens_60/sTime_30/sampleTime_1.000000e-01
%
%   The function recursively searches baseDir for output_database_1.db
%   and selects the newest match.
%
% INPUTS
%   baseDir  (char | string, optional) — path to the WiLabVISim run folder
%                                        containing the SQLite database.
%
% OUTPUT
%   wilab_export.mat  — saved in the same folder as the database; contains:
%     dbfile          (char)   — full path to the source database
%     tbls            (string) — list of tables found in the database
%     vLoc_table      (table)  — vehicle location table (vLocation)
%     vLoc_mat        (double) — vLoc_table as a numeric matrix
%     vLoc_cols       (string) — column names for vLoc_mat
%     intf_table      (table)  — interference path table (firstReflectLog)
%     intf_mat        (double) — intf_table as a numeric matrix
%     intf_cols       (string) — column names for intf_mat
%     runInfo_table   (table)  — simulation metadata (if present)
%
% REQUIRED WILAB SCHEMA
% ---------------------
%   vLocation table columns:
%     time, Vid, V_x, V_y, V_speed,
%     corner_p1_x, corner_p1_y, corner_p2_x, corner_p2_y,
%     corner_p3_x, corner_p3_y, corner_p4_x, corner_p4_y
%
%   firstReflectLog table columns:
%     time, interferedId, interferedSensor_x, interferedSensor_y,
%     interferingId, interferingSensor_x, interferingSensor_y,
%     lengthPathFirst, lengthPathSecond, reflectPointX, reflectPointY
%
% DEPENDENCIES
%   Requires MATLAB Database Toolbox (sqlite function).
%
% AUTHORS
%   Arianit Preniqi, Oliver Lang, Stefan Schmalzl, Reinhard Feger
%
%   Johannes Kepler University Linz, Linz, Austria
% VERSION
%   v1.0  —  26.03.2026
%
% SEE ALSO
%   main_sim, main_sim_FF
 
% =========================================================================
 
baseDir = fullfile(pwd,'Output\antenna4\dens_60\sTime_30\sampleTime_1.000000e-01');
 
% --- find newest output_database_*.db (works in run folder or any root under Output)
dbList = dir(fullfile(baseDir, '**', 'output_database_1.db'));
assert(~isempty(dbList), 'No output_database_*.db found under: %s', baseDir);
[~, idx] = max([dbList.datenum]);
dbfile = fullfile(dbList(idx).folder, dbList(idx).name);
fprintf('DB: %s\n', dbfile);
 
% --- connect
conn = sqlite(dbfile,'readonly');
 
% --- list tables
tblsQ = fetch(conn, "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;");
tbls  = toStringVector(tblsQ, 'name');
fprintf('Tables: %s\n', strjoin(cellstr(tbls), ', '));
 
% --- exact table names used by WiLabVISim
vLocName = firstMatch(tbls, ["vLocation"]);
intfName = firstMatch(tbls, ["firstReflectLog"]);
runInfo  = firstMatch(tbls, ["runInfo"]);
 
% --- read tables (empty table if missing)
vLoc_table    = readSqlTable(conn, vLocName);
intf_table    = readSqlTable(conn, intfName);
runInfo_table = readSqlTable(conn, runInfo);
 
% quick stats
if ~isempty(vLocName)
    nV = fetchScalar(conn, sprintf('SELECT COUNT(*) FROM "%s";', vLocName));
    tV = fetch(conn, sprintf('SELECT MIN(time), MAX(time) FROM "%s";', vLocName));
    [tVmin,tVmax] = twoScalars(tV);
    fprintf('vLocation: %d rows | t=[%.3f .. %.3f]\n', nV, tVmin, tVmax);
else
    fprintf('vLocation table not present.\n');
end
if ~isempty(intfName)
    nI = fetchScalar(conn, sprintf('SELECT COUNT(*) FROM "%s";', intfName));
    tI = fetch(conn, sprintf('SELECT MIN(time), MAX(time) FROM "%s";', intfName));
    [tImin,tImax] = twoScalars(tI);
    fprintf('firstReflectLog: %d rows | t=[%.3f .. %.3f]\n', nI, tImin, tImax);
else
    fprintf('firstReflectLog table not present.\n');
end
if ~isempty(runInfo)
    fprintf('runInfo present (sim flags/metadata).\n');
end
 
close(conn);
 
% --- convert to numeric matrices using WiLabVISim schema
[vLoc_mat, vLoc_cols] = vloc_to_matrix_wilab(vLoc_table);
[intf_mat, intf_cols] = intf_to_matrix_wilab(intf_table);
 
% --- save
outMat = fullfile(fileparts(dbfile), 'wilab_export.mat');
save(outMat, 'dbfile','tbls','vLoc_table','vLoc_mat','vLoc_cols','intf_table','intf_mat','intf_cols','runInfo_table','-v7.3');
fprintf('Saved: %s\n', outMat);
 
% --- peek
if ~isempty(vLoc_mat)
    fprintf('vLoc_mat: %d x %d | cols: %s\n', size(vLoc_mat,1), size(vLoc_mat,2), strjoin(cellstr(vLoc_cols), ', '));
else
    fprintf('vLoc_mat is empty (did logging fire? sampleTime < simulationTime? Correct table names?).\n');
end
if ~isempty(intf_mat)
    fprintf('intf_mat: %d x %d | cols: %s\n', size(intf_mat,1), size(intf_mat,2), strjoin(cellstr(intf_cols), ', '));
else
    fprintf('intf_mat is empty (did logPotentialInterferers run?).\n');
end
end
 
% =========================================================================
%% ── Local Helper Functions ───────────────────────────────────────────────
% =========================================================================
 
function s = toStringVector(x, preferCol)
    if istable(x)
        if nargin>=2 && ismember(preferCol, x.Properties.VariableNames)
            s = string(x.(preferCol));
        else
            vn = x.Properties.VariableNames{1};
            s = string(x.(vn));
        end
    elseif iscell(x)
        s = string(x(:));
    elseif isstring(x)
        s = x(:);
    else
        s = string(x(:));
    end
end
 
function name = firstMatch(names, candidates)
    name = "";
    if isempty(names), return; end
    L = lower(string(names));
    for c = string(candidates)
        idx = find(L == lower(c), 1);
        if ~isempty(idx), name = names(idx); return; end
    end
    for c = string(candidates)
        idx = find(contains(L, lower(c)), 1);
        if ~isempty(idx), name = names(idx); return; end
    end
end
 
function T = readSqlTable(conn, tname)
    if tname == "" || isempty(tname), T = table(); return; end
    rows = fetch(conn, sprintf('SELECT * FROM "%s";', tname));
    if istable(rows), T = rows; return; end
    pragma = fetch(conn, sprintf('PRAGMA table_info("%s");', tname));
    colnames = extractColnames(pragma);
    if isempty(rows)
        if isempty(colnames), T = table();
        else, T = cell2table(cell(0,numel(colnames)),'VariableNames',cellstr(colnames));
        end
    else
        if isempty(colnames), colnames = "Var"+(1:size(rows,2)); end
        T = cell2table(rows,'VariableNames',cellstr(colnames));
    end
end
 
function names = extractColnames(pragma)
    names = strings(0,1);
    if istable(pragma)
        if any(strcmpi(pragma.Properties.VariableNames,'name')), names = string(pragma.('name'));
        elseif width(pragma) >= 2, names = string(pragma{:,2}); end
    elseif iscell(pragma) && ~isempty(pragma)
        names = string(pragma(:,2));
    end
end
 
function x = fetchScalar(conn, q)
    r = fetch(conn,q);
    if iscell(r), x = r{1};
    elseif istable(r), x = r{1,1};
    else, x = r(1);
    end
end
 
function [a,b] = twoScalars(r)
    if iscell(r), a=r{1}; b=r{2};
    elseif istable(r), a=r{1,1}; b=r{1,2};
    else, a=r(1); b=r(2);
    end
    if isempty(a) || ~isfinite(a), a = NaN; end
    if isempty(b) || ~isfinite(b), b = NaN; end
end
 
function [M, cols] = vloc_to_matrix_wilab(T)
    M = []; cols = strings(0,1);
    if isempty(T), return; end
    must = ["time","Vid","V_x","V_y","V_speed", ...
            "corner_p1_x","corner_p1_y","corner_p2_x","corner_p2_y", ...
            "corner_p3_x","corner_p3_y","corner_p4_x","corner_p4_y"];
    if ~all(ismember(must, string(T.Properties.VariableNames)))
        warning('vLocation does not match expected WiLab schema. Found: %s', strjoin(T.Properties.VariableNames, ', '));
        return;
    end
    M = [double(T.time), double(T.Vid), double(T.V_x), double(T.V_y), double(T.V_speed), ...
         double(T.corner_p1_x), double(T.corner_p1_y), double(T.corner_p2_x), double(T.corner_p2_y), ...
         double(T.corner_p3_x), double(T.corner_p3_y), double(T.corner_p4_x), double(T.corner_p4_y)];
    cols = ["time","Vid","V_x","V_y","V_speed", ...
            "corner_p1_x","corner_p1_y","corner_p2_x","corner_p2_y", ...
            "corner_p3_x","corner_p3_y","corner_p4_x","corner_p4_y"];
end
 
function [M, cols] = intf_to_matrix_wilab(R)
    M = []; cols = strings(0,1);
    if isempty(R), return; end
    need = ["time","interferedId","interferedSensor_x","interferedSensor_y", ...
            "interferingId","interferingSensor_x","interferingSensor_y", ...
            "lengthPathFirst","lengthPathSecond","reflectPointX","reflectPointY"];
    if ~all(ismember(need, string(R.Properties.VariableNames)))
        warning('firstReflectLog does not match expected WiLab schema. Found: %s', strjoin(R.Properties.VariableNames, ', '));
        return;
    end
    M = [double(R.time), double(R.interferedId), double(R.interferingId), ...
         double(R.interferedSensor_x), double(R.interferedSensor_y), ...
         double(R.interferingSensor_x), double(R.interferingSensor_y), ...
         double(R.lengthPathFirst), double(R.lengthPathSecond), ...
         double(R.reflectPointX), double(R.reflectPointY)];
    cols = ["time","victim_id","source_id","rx_x","rx_y","tx_x","tx_y","L1","L2","rpX","rpY"];
end
