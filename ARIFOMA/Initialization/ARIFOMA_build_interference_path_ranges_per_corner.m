function pathByCorner = ARIFOMA_build_interference_path_ranges_per_corner(links, vehId, corners)
% ARIFOMA_BUILD_INTERFERENCE_PATH_RANGES_PER_CORNER
%   Extract LOS and reflected interference path distance arrays per victim
%   corner from the ARIFOMA interference link table.
%
% PURPOSE
% -------
%   Converts the flat interference link table (one row per path) into a
%   per-corner struct that groups path distances by type (LOS / reflected)
%   and corner label.  The resulting pathByCorner struct is consumed by
%   ARIFOMA_build_engine_inputs to populate engine.dist for each corner's
%   CPI processing step.
%
%   Only paths where VictimVehID matches vehId are retained.  All other
%   rows are discarded.
%
% INPUTS
%   links    (table)         — interference link table produced by
%                              ARIFOMA_build_interference_link_table_from_frame;
%                              required columns:
%                                VictimVehID  (double)
%                                VictimCorner (string)
%                                PathType     (string) — "LOS" or "REFL"
%                                d1_m         (double) — first leg length [m]
%                                d2_m         (double) — second leg length [m]
%   vehId    (double)        — victim vehicle ID; only links with
%                              VictimVehID == vehId are processed
%   corners  (string array)  — corner labels to build entries for,
%                              e.g. ["FL","FR","RL","RR"]
%
% OUTPUT
%   pathByCorner  (struct) — struct indexed by corner name (e.g. pathByCorner.FL)
%                            each entry is a struct with fields:
%                              LOS       (double, 1×N_LOS)  — LOS path distances [m]
%                              REFTotal  (double, 1×N_REF)  — total reflected path
%                                                             lengths (d1+d2) [m]
%                              d1        (double, 1×N_REF)  — first reflected leg [m]
%                              d2        (double, 1×N_REF)  — second reflected leg [m]
%                            All arrays are row vectors.  Arrays are empty
%                            when no paths of that type exist for that corner.
%
% NOTES
% -----
%   - Non-finite or negative distance values are removed before storage to
%     prevent downstream INR / power calculations from receiving bad inputs.
%   - For reflected paths, d1 and d2 arrays are truncated to the same
%     length before computing REFTotal = d1 + d2 to guard against
%     mismatched row counts from upstream table construction.
%   - If links is empty or not a table, all corners are initialised with
%     empty arrays and the function returns immediately.
%
% SEE ALSO
%   ARIFOMA_build_interference_link_table_from_frame,
%   ARIFOMA_build_engine_inputs,
%   ARIFOMA_prepare_scene_views

% =========================================================================
    pathByCorner = struct();
if isempty(links) || ~istable(links)
for c = 1:numel(corners)
            cn = char(corners(c));
            pathByCorner.(cn) = struct('LOS',[],'REFTotal',[],'d1',[],'d2',[]);
end
return
end
% keep only this victim vehicle
    L = links(links.VictimVehID == vehId, :);
for c = 1:numel(corners)
        cn = char(corners(c));
        Lc = L(L.VictimCorner == cn, :);
% LOS rows: use d1_m as the straight range
        los = double(Lc.d1_m(Lc.PathType == "LOS"));
        los = los(isfinite(los) & los>=0);
% Reflected rows: total path length d1+d2
        d1  = double(Lc.d1_m(Lc.PathType == "REFL"));
        d2  = double(Lc.d2_m(Lc.PathType == "REFL"));
        d1  = d1(isfinite(d1) & d1>=0);
        d2  = d2(isfinite(d2) & d2>=0);
        n   = min(numel(d1), numel(d2));
        reflTot = d1(1:n) + d2(1:n);
        pathByCorner.(cn) = struct( ...
'LOS',      los(:).', ...
'REFTotal', reflTot(:).', ...
'd1',       d1(:).', ...
'd2',       d2(:).' );
end
end