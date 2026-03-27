function links = ARIFOMA_build_interference_link_table_from_frame( ...
        frame, vLoc_table, Param, tol, cfg)
% ARIFOMA_BUILD_INTERFERENCE_LINK_TABLE_FROM_FRAME
%   Build a per-path interference link table from a pre-extracted scene frame.
%
% PURPOSE
% -------
%   Converts the raw WiLabVISim interference path records stored in
%   frame.R_draw into a structured MATLAB table where each row describes
%   one directed interference link between a victim sensor corner and an
%   interfering sensor corner, annotated with path type and distances.
%
%   The function resolves physical sensor corner positions (from vLoc_table
%   corner geometry) to logical corner labels (FL / FR / RR / RL) and
%   looks up the corresponding SensorID from the ARIFOMA parameter table.
%
% INPUTS
%   frame       (struct)  — scene frame struct produced by
%                           ARIFOMA_frame_extract_mc; required fields:
%                             frame.R_draw  — interference path table
%                                             (WiLabVISim firstReflectLog)
%                             frame.t0      — snapshot timestamp [s]
%   vLoc_table  (table)   — vehicle location table with columns:
%                             time, Vid,
%                             corner_p1_x/y … corner_p4_x/y
%   Param       (table)   — ARIFOMA parameter table produced by
%                           ARIFOMA_build_parameter_table; must contain
%                           columns VehID, Corner, and optionally SensorID
%   tol         (double)  — time tolerance for slicing vLoc_table at
%                           frame.t0 [s]  (default: 0.05)
%   cfg         (struct)  — global configuration struct; recognised field:
%                             cfg.FrontSensor (logical, default false)
%                               true  → single FF sensor layout
%                               false → four corner sensors FL/FR/RR/RL
%
% OUTPUT
%   links  (table) — interference link table with one row per path in
%                    frame.R_draw that could be fully resolved.
%                    Columns:
%                      VictimVehID        (double) — victim vehicle ID
%                      VictimCorner       (string) — victim corner label
%                      VictimSensorID     (double) — victim sensor ID (NaN if absent)
%                      InterfererVehID    (double) — interferer vehicle ID
%                      InterfererCorner   (string) — interferer corner label
%                      InterfererSensorID (double) — interferer sensor ID (NaN if absent)
%                      PathType           (string) — "LOS" or "REFL"
%                      d1_m               (double) — first path leg length [m]
%                      d2_m               (double) — second path leg length [m];
%                                                    0 for LOS paths
%
% CORNER LABEL MAPPING
% --------------------
%   Four-corner mode (cfg.FrontSensor = false, default):
%     Corner labels are resolved by finding the vehicle corner position
%     (corner_p1 … corner_p4) nearest to the sensor position reported in
%     frame.R_draw, then mapping the index to a label:
%       Index 1 → FL  |  Index 2 → FR  |  Index 3 → RR  |  Index 4 → RL
%
%   Front-sensor mode (cfg.FrontSensor = true):
%     Only one sensor exists per vehicle. All paths are assigned the
%     label "FF" directly — no nearest-corner disambiguation needed.
%
% PATH TYPE CLASSIFICATION
% ------------------------
%   A path is classified as LOS if any of the following hold:
%     - lengthPathSecond is non-finite or zero
%     - reflectPointX or reflectPointY is non-finite
%   Otherwise the path is classified as REFL (first-order reflected).
%
% NOTES
% -----
%   - Rows in frame.R_draw for which either the victim or interferer vehicle
%     cannot be found in vLoc_table at time frame.t0 are silently skipped.
%   - SensorID lookup returns NaN if the SensorID column is absent from
%     Param or if no matching (VehID, Corner) row exists.
%   - The output table is empty (zero rows) if frame.R_draw is absent,
%     empty, or all rows are skipped.
%
% SEE ALSO
%   ARIFOMA_frame_extract_mc,
%   ARIFOMA_prepare_scene_views,
%   ARIFOMA_build_parameter_table
 
% =========================================================================
 
%% ── Default Arguments ────────────────────────────────────────────────────
if nargin < 4,  tol = 0.05;  end   % [s] time tolerance for vLoc slice
if nargin < 5,  cfg = struct(); end
 
% Determine sensor layout mode from cfg.FrontSensor:
%   true  → single front-facing sensor, corner label "FF", one position point
%   false → four corner sensors FL/FR/RR/RL (default)
if isfield(cfg, 'FrontSensor') && cfg.FrontSensor
    CornerMap  = ["FF"];
    numCorners = 1;
else
    CornerMap  = ["FL", "FR", "RR", "RL"];
    numCorners = 4;
end
 
%% ── Pre-allocate Empty Output Table ─────────────────────────────────────
%   Pre-allocation avoids schema inference from the first appended row and
%   guarantees the correct column types even when the table is returned empty.
 
VictimVehID        = double.empty(0,1);
VictimCorner       = strings(0,1);
VictimSensorID     = double.empty(0,1);
InterfererVehID    = double.empty(0,1);
InterfererCorner   = strings(0,1);
InterfererSensorID = double.empty(0,1);
PathType           = strings(0,1);
d1_m               = double.empty(0,1);
d2_m               = double.empty(0,1);
 
links = table( ...
    VictimVehID, VictimCorner, VictimSensorID, ...
    InterfererVehID, InterfererCorner, InterfererSensorID, ...
    PathType, d1_m, d2_m);
 
%% ── Early Return if No Paths Available ──────────────────────────────────
if ~isfield(frame, 'R_draw') || isempty(frame.R_draw)
    return;
end
 
%% ── Slice Vehicle Locations at Snapshot Time ─────────────────────────────
%   Retain only vehicle rows whose timestamp is within tol of frame.t0.
%   getVrow is a convenience handle for looking up a single vehicle row.
 
R      = frame.R_draw;
rowsV  = abs(double(vLoc_table.time) - frame.t0) <= tol;
Vt     = vLoc_table(rowsV, :);
 
%% ── Corner Label Map ─────────────────────────────────────────────────────
%   CornerMap and numCorners are set above in the Default Arguments section
%   based on cfg.FrontSensor. No hardcoding here.
%     Four-corner mode: 1→FL | 2→FR | 3→RR | 4→RL
%     Front-sensor mode: always "FF"
 
%% ── Lookup Handles ───────────────────────────────────────────────────────
getVrow    = @(vid) Vt(find(double(Vt.Vid) == double(vid), 1, 'first'), :);
corner_idx = @(vr, sx, sy) local_corner_idx(vr, sx, sy, numCorners);
 
%% ── Main Path Loop ───────────────────────────────────────────────────────
%   Each row in R corresponds to one interference path (LOS or reflected).
%   For each row:
%     1. Resolve the victim vehicle row from vLoc_table.
%     2. Map the victim sensor position to a corner label.
%     3. Look up the victim SensorID from Param.
%     4. Repeat steps 1-3 for the interferer.
%     5. Classify the path as LOS or REFL.
%     6. Append the link to the output table.
 
for i = 1:height(R)
 
    % ── Victim Endpoint Resolution ────────────────────────────────────────
    vID  = R.interferedId(i);
    vrx  = getVrow(vID);
    if isempty(vrx),  continue;  end   % vehicle not in vLoc slice — skip
 
    k_v  = corner_idx(vrx, R.interferedSensor_x(i), R.interferedSensor_y(i));
    labV = CornerMap(k_v);
 
    % Look up victim SensorID from Param (NaN if column absent or no match).
    sidV = NaN;
    if ismember('SensorID', Param.Properties.VariableNames)
        rowV = Param(Param.VehID == vID & Param.Corner == labV, :);
        if ~isempty(rowV),  sidV = rowV.SensorID(1);  end
    end
 
    % ── Interferer Endpoint Resolution ───────────────────────────────────
    iID  = R.interferingId(i);
    vtx  = getVrow(iID);
    if isempty(vtx),  continue;  end   % vehicle not in vLoc slice — skip
 
    k_i  = corner_idx(vtx, R.interferingSensor_x(i), R.interferingSensor_y(i));
    labI = CornerMap(k_i);
 
    % Look up interferer SensorID from Param.
    sidI = NaN;
    if ismember('SensorID', Param.Properties.VariableNames)
        rowI = Param(Param.VehID == iID & Param.Corner == labI, :);
        if ~isempty(rowI),  sidI = rowI.SensorID(1);  end
    end
 
    % ── Path Type Classification ──────────────────────────────────────────
    %   A path is LOS if the second leg length is absent/zero or the
    %   reflection point coordinates are non-finite.
    L2    = R.lengthPathSecond(i);
    isLOS = ~isfinite(L2) || L2 <= 0 || ...
            ~isfinite(R.reflectPointX(i)) || ~isfinite(R.reflectPointY(i));
 
    ptype = string(ifelse(isLOS, 'LOS', 'REFL'));
    d1    = R.lengthPathFirst(i);
    d2    = ifelse(isLOS, 0, max(0, L2));   % second leg is 0 for LOS
 
    % ── Append Link Row ───────────────────────────────────────────────────
    links = [links; {vID, labV, sidV, iID, labI, sidI, ptype, d1, d2}];   %#ok<AGROW>
 
end   % path loop
 
end   % ARIFOMA_build_interference_link_table_from_frame
 
% =========================================================================
%% ── Local Helper Functions ───────────────────────────────────────────────
% =========================================================================
 
% -------------------------------------------------------------------------
function k = local_corner_idx(vr, sx, sy, numCorners)
% LOCAL_CORNER_IDX
%   Find the vehicle corner nearest to a given sensor position.
%
%   Four-corner mode (numCorners = 4):
%     Computes the Euclidean distance from the sensor position (sx, sy) to
%     each of the four vehicle corner positions (corner_p1 … corner_p4) and
%     returns the index of the nearest corner, clamped to [1, 4].
%       maps to: 1→FL, 2→FR, 3→RR, 4→RL
%
%   Front-sensor mode (numCorners = 1):
%     Only one sensor exists per vehicle — always returns 1 (→ "FF").
%     No geometric lookup needed.
%
%   INPUTS
%     vr         (table row) — single vehicle row from vLoc_table
%     sx         (double)    — sensor x-coordinate [m]
%     sy         (double)    — sensor y-coordinate [m]
%     numCorners (int)       — 1 for FF mode, 4 for four-corner mode
%
%   OUTPUT
%     k  (int) — nearest corner index; always 1 in FF mode
 
    % Front-sensor mode: only one possible label, skip geometry lookup.
    if numCorners == 1
        k = 1;
        return;
    end
 
    % Four-corner mode: find nearest corner by Euclidean distance.
    C = [vr.corner_p1_x, vr.corner_p1_y; ...
         vr.corner_p2_x, vr.corner_p2_y; ...
         vr.corner_p3_x, vr.corner_p3_y; ...
         vr.corner_p4_x, vr.corner_p4_y];
 
    [~, k] = min(hypot(double(C(:,1)) - double(sx), ...
                       double(C(:,2)) - double(sy)));
 
    % Clamp to valid range as a safety guard.
    k = max(1, min(4, k));
end
 
% -------------------------------------------------------------------------
function y = ifelse(cond, a, b)
% IFELSE  Inline conditional expression (ternary operator equivalent).
%
%   Returns a if cond is true, b otherwise.  Used to keep path-type and
%   distance assignment on a single readable line.
%
%   INPUTS
%     cond  (logical) — condition to evaluate
%     a               — value returned when cond is true
%     b               — value returned when cond is false
%
%   OUTPUT
%     y               — a or b
 
    if cond,  y = a;  else,  y = b;  end
end