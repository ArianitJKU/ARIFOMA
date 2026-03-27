function ARIFOMA_draw_scene_snapshot(vLoc, intf, t_query, tol, maxPaths, stretchY, axTop, varargin)
% ARIFOMA_DRAW_SCENE_SNAPSHOT
%   Render one scenario snapshot into a MATLAB axes handle.
%
% PURPOSE
% -------
%   This function is a pure visualizer — it does not modify any simulation
%   state.  Given a vehicle-location table and an interference-path table,
%   it slices both datasets at the nearest available timestamp to t_query
%   and produces a top-down highway scene showing:
%
%     - All vehicles rendered as speed-colored car icons or polygons
%     - LOS interference paths (solid red lines)
%     - First-order reflected interference paths (dashed grey lines)
%     - Reflection points (yellow markers)
%     - Optional focus highlight (blue circle) around a selected vehicle
%     - Smooth camera panning that follows a focus vehicle across frames
%
%   The function is designed to be called repeatedly inside a frame loop
%   (e.g. a video export loop) without accumulating graphics objects or
%   triggering MATLAB colorbar listener crashes.
%
% VIDEO-SAFE DESIGN NOTES
% -----------------------
%   - The colorbar and colormap are created ONCE per axis on the first call
%     and reused on all subsequent calls via appdata caching.
%   - axis(ax,'manual') is never called inside the frame loop as this
%     triggers a MATLAB listener crash when a colorbar is present.
%   - Axis limits are set exclusively via set(ax,'XLim',...,'YLim',...).
%   - Camera panning uses exponential smoothing (panAlpha) to avoid
%     abrupt jumps between frames.
%
% INPUTS (positional)
%   vLoc      — vehicle location data; either:
%                 (a) a table with columns: time, Vid, V_x, V_y, V_speed,
%                     corner_p1_x/y … corner_p4_x/y  (WiLabVISim schema), or
%                 (b) a numeric matrix with ≥13 columns:
%                     [t, Vid, Vx, Vy, speed, p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y]
%   intf      — interference path table (WiLabVISim firstReflectLog schema)
%   t_query   — requested snapshot time [s]; nearest available timestamp
%               is selected automatically
%   tol       — time tolerance around the selected t0 [s]  (default: 0.2)
%   maxPaths  — maximum number of interference paths to render (default: 200)
%   stretchY  — vertical display stretch factor, clamped to [0.1, 1]
%               (default: 0.35)
%   axTop     — target MATLAB axes handle
%
% NAME-VALUE OPTIONS
%   All options use manual parsing (no inputParser) for speed.
%
%   'VehicleId'        [] or numeric — vehicle ID to focus / follow;
%                                      [] disables focus mode
%   'Role'             string        — role filter for path display:
%                                      'victim' | 'source' | 'either'
%                                      (default: 'either')
%   'ShowAllVehicles'  logical       — render all vehicles when true;
%                                      when false only vehicles involved
%                                      in interference paths are shown
%                                      (default: true)
%   'WinX'             double        — camera window half-width [m]
%                                      (default: 400)
%   'PanAlpha'         double ∈[0,1] — exponential smoothing factor for
%                                      camera panning; 0 = frozen,
%                                      1 = instant snap  (default: 0.25)
%   'UseCarIcons'      logical       — render chamfered car glyphs (true)
%                                      or raw corner polygons (false)
%                                      (default: true)
%   'NumLanes'         int           — number of lanes; used for lane
%                                      direction assignment (default: 6)
%   'LanesWest'        int array     — lane indices assigned westbound
%                                      heading (default: [1 2 3])
%   'CarLenClamp'      [min max]     — vehicle length clamp range [m]
%                                      (default: [3.6 5.0])
%   'CarWidClamp'      [min max]     — vehicle width clamp range [m]
%                                      (default: [1.6 2.2])
%   'DebugCarHeading'  logical       — overlay quiver heading arrows on
%                                      the first 8 vehicles (default: false)
%
% OUTPUTS
%   None — all output is rendered directly into axTop.
%
% NOTES
%   - Speed colorbar is hard-clamped to [100 160] km/h to match the
%     highway scenario speed range used in the published case study.
%     Adjust the clim() call to change this range.
%   - Lane centres are estimated via k-means on the first call and cached
%     in axTop appdata for stability across frames.
%   - Interference paths are downselected to maxPaths by total path length
%     (longest paths rendered first) when the table exceeds maxPaths rows.
%
% SEE ALSO
%   ARIFOMA_frame_extract_mc, ARIFOMA_frame_draw_mc, main_sim

% =========================================================================

%% ── Positional Argument Defaults ─────────────────────────────────────────
if nargin < 4 || isempty(tol),      tol      = 0.2;   end
if nargin < 5 || isempty(maxPaths), maxPaths = 200;   end
if nargin < 6 || isempty(stretchY), stretchY = 0.35;  end
stretchY = max(0.1, min(1, stretchY));   % clamp to valid range

%% ── Name-Value Option Defaults ───────────────────────────────────────────
vehIdOpt    = [];
roleOpt     = "either";
showAllVeh  = true;
WinX        = 400;
panAlpha    = 0.25;
useCars     = true;
NumLanes    = 6;
LanesWest   = [1 2 3];
CarLenClamp = [3.6 5.0];
CarWidClamp = [1.6 2.2];
dbgHead     = false;

%% ── Name-Value Manual Parsing ────────────────────────────────────────────
%   inputParser is deliberately avoided here for performance inside a
%   frame loop.  Options must be supplied as ('Name', value) pairs.

if mod(numel(varargin), 2) ~= 0
    error('ARIFOMA:drawSceneSnapshot:oddNameValue', ...
          'Name-value options must be supplied in pairs.');
end

for k = 1:2:numel(varargin)
    name = lower(string(varargin{k}));
    val  = varargin{k+1};

    switch name
        case "vehicleid"
            vehIdOpt = val;
        case "role"
            roleOpt = lower(string(val));
        case "showallvehicles"
            showAllVeh = logical(val);
        case "winx"
            WinX = max(1, double(val));
        case "panalpha"
            panAlpha = max(0, min(1, double(val)));
        case "usecaricons"
            useCars = logical(val);
        case "numlanes"
            NumLanes = max(2, round(val));
        case "laneswest"
            LanesWest = unique(val(:))';
        case "carlenclamp"
            CarLenClamp = double(val(:))';
        case "carwidclamp"
            CarWidClamp = double(val(:))';
        case "debugcarheading"
            dbgHead = logical(val);
        otherwise
            error('ARIFOMA:drawSceneSnapshot:unknownOption', ...
                  'Unknown option: "%s".', name);
    end
end

% Validate role string.
if ~ismember(roleOpt, ["victim","source","either"])
    error('ARIFOMA:drawSceneSnapshot:invalidRole', ...
          'Role must be "victim", "source", or "either".');
end

%% ── Input Validation ─────────────────────────────────────────────────────
isVTable = istable(vLoc);

if ~isVTable
    if ~isnumeric(vLoc) || size(vLoc, 2) < 13
        error('ARIFOMA:drawSceneSnapshot:invalidVLoc', ...
              'Numeric vLoc must have at least 13 columns.');
    end
end

if ~istable(intf)
    error('ARIFOMA:drawSceneSnapshot:invalidIntf', ...
          'intf must be a table.');
end

% Normalise interference table column names to a canonical set.
intf = normalize_intf_columns(intf);

%% ── Timestamp Selection ──────────────────────────────────────────────────
%   The nearest available timestamp across both tables is selected.
%   This avoids errors when vLoc and intf are sampled at slightly
%   different times.

if isVTable
    tv = unique(double(vLoc.time));
else
    tv = unique(double(vLoc(:,1)));
end

ti   = unique(double(intf.time));
allT = unique([tv(:); ti(:)]);

if isempty(allT)
    error('ARIFOMA:drawSceneSnapshot:noTimestamps', ...
          'No timestamps found in vLoc or intf.');
end

[~, kmin] = min(abs(allT - t_query));
t0        = allT(kmin);

%% ── Data Slicing at t0 ───────────────────────────────────────────────────
%   Select rows within ±tol of t0 from both tables.

if isVTable
    rowsV = abs(double(vLoc.time) - t0) <= tol;
    V = vLoc(rowsV, :);
    M = [];
else
    rowsV = abs(double(vLoc(:,1)) - t0) <= tol;
    M = vLoc(rowsV, :);
    V = table();
end

rowsI = abs(double(intf.time) - t0) <= tol;
R     = intf(rowsI, :);

if (isVTable && isempty(V)) || (~isVTable && isempty(M))
    error('ARIFOMA:drawSceneSnapshot:noVehicles', ...
          'No vehicle rows found near t=%.3f s (tol=%.3f s).', t0, tol);
end

%% ── Interference Row Sanitisation ───────────────────────────────────────
%   Negative placeholder values in WiLabVISim logs indicate missing data.
%   These are replaced with NaN and the corresponding rows are removed.

if ~isempty(R)
    colsToSanitize = intersect( ...
        {'interferedSensor_x','interferedSensor_y', ...
         'interferingSensor_x','interferingSensor_y', ...
         'reflectPointX','reflectPointY','lengthPathSecond'}, ...
        R.Properties.VariableNames);

    for c = 1:numel(colsToSanitize)
        v = double(R.(colsToSanitize{c}));
        v(v <= -0.5) = NaN;
        R.(colsToSanitize{c}) = v;
    end

    % Remove rows where either endpoint coordinate is undefined.
    bad = isnan(R.interferedSensor_x)  | isnan(R.interferedSensor_y) | ...
          isnan(R.interferingSensor_x) | isnan(R.interferingSensor_y);
    R(bad, :) = [];
end

%% ── Focus Filter (Interference Paths) ───────────────────────────────────
%   When a focus vehicle is specified, retain only paths where that vehicle
%   appears in the role defined by roleOpt.

if ~isempty(vehIdOpt)
    switch roleOpt
        case "victim"
            R = R(R.interferedId  == vehIdOpt, :);
        case "source"
            R = R(R.interferingId == vehIdOpt, :);
        otherwise
            R = R(R.interferedId  == vehIdOpt | ...
                  R.interferingId == vehIdOpt, :);
    end
end

%% ── Unified Vehicle Table Construction ───────────────────────────────────
%   Both vLoc_table and numeric vLoc_mat are normalised into a single
%   table S with a consistent column schema for downstream rendering.

if isVTable
    need = ["Vid","V_x","V_y","V_speed", ...
            "corner_p1_x","corner_p1_y","corner_p2_x","corner_p2_y", ...
            "corner_p3_x","corner_p3_y","corner_p4_x","corner_p4_y"];
    S          = V(:, need);
    S.V_speed  = S.V_speed * 3.6;   % convert m/s → km/h for colorbar
else
    S             = table();
    S.Vid         = M(:,2);
    S.V_x         = M(:,3);
    S.V_y         = M(:,4);
    S.V_speed     = M(:,5);
    S.corner_p1_x = M(:,6);  S.corner_p1_y = M(:,7);
    S.corner_p2_x = M(:,8);  S.corner_p2_y = M(:,9);
    S.corner_p3_x = M(:,10); S.corner_p3_y = M(:,11);
    S.corner_p4_x = M(:,12); S.corner_p4_y = M(:,13);
end

% When ShowAllVehicles = false, retain only vehicles involved in paths.
if ~isempty(vehIdOpt) && ~showAllVeh
    involved = ~isempty(R) && unique([R.interferedId; R.interferingId]) || vehIdOpt;
    S        = S(ismember(S.Vid, involved), :);
end

%% ── Axes Initialisation (Colorbar Caching) ───────────────────────────────
%   The colorbar and colormap are created ONCE on the first call and cached
%   in the axes appdata.  Recreating them every frame causes instability
%   with MATLAB's colorbar listener mechanism.

ax = axTop;

if ~isappdata(ax, 'viz_initialized')
    setappdata(ax, 'viz_initialized', true);
    setappdata(ax, 'viz_cmap', parula(256));

    hcb = colorbar(ax, 'Location','eastoutside');
    setappdata(ax, 'viz_colorbar', hcb);

    % Set limit modes once to prevent automatic rescaling on cla().
    set(ax, 'XLimMode','manual', 'YLimMode','manual');
end

cmap = getappdata(ax, 'viz_cmap');
colormap(ax, cmap);
hcb = getappdata(ax, 'viz_colorbar');

% Clear axes content without destroying the cached colorbar.
cla(ax);
grid(ax, 'off');
box(ax,  'on');
set(ax, 'NextPlot','add', 'Layer','top');

%% ── Legend Proxy Objects ─────────────────────────────────────────────────
%   Invisible proxy objects are created after cla() to populate the legend
%   without depending on actual path or vehicle graphics objects.

hVeh = patch(ax, NaN, NaN, 'w', ...
    'HandleVisibility','on', 'DisplayName','Vehicles');
hLoS = plot(ax, NaN, NaN, '-', 'Color',[0.85 0 0], 'LineWidth',2, ...
    'HandleVisibility','on', 'DisplayName','LOS Interference Path');
hRef = plot(ax, NaN, NaN, '--', 'Color',[0.5 0.5 0.5], 'LineWidth',2, ...
    'HandleVisibility','on', 'DisplayName','First-order Reflected Path');
hRP  = plot(ax, NaN, NaN, 'o', 'MarkerSize',6, ...
    'MarkerFaceColor','y', 'MarkerEdgeColor','k', ...
    'HandleVisibility','on', 'DisplayName','Reflection Point');

if isempty(vehIdOpt)
    legend(ax, [hVeh hLoS hRef hRP], ...
        'Location','eastoutside', 'Color','w', 'Box','on');
else
    hFocusProxy = plot(ax, NaN, NaN, 'b-', 'LineWidth',4, ...
        'HandleVisibility','on', 'DisplayName','Victim Vehicle');
    legend(ax, [hVeh hLoS hRef hRP hFocusProxy], ...
        'Location','northoutside', 'Color','w', 'Box','on');
end

%% ── Camera Panning ───────────────────────────────────────────────────────
%   The camera follows the focus vehicle using exponential smoothing.
%   The previous camera centre is cached in appdata between frames.

halfW = WinX / 2;

xc_target = [];
if ~isempty(vehIdOpt)
    rFocus = S(S.Vid == vehIdOpt, :);
    if ~isempty(rFocus)
        xc_target = double(rFocus.V_x(1));
    end
end
if isempty(xc_target)
    xc_target = median(double(S.V_x));
end

% Retrieve or initialise the previous camera centre.
if isappdata(ax, 'viz_last_xc')
    xc_prev = getappdata(ax, 'viz_last_xc');
else
    xc_prev = xc_target;
end

% Exponential smoothing: xc = (1-α)*xc_prev + α*xc_target
xc = (1 - panAlpha) * xc_prev + panAlpha * xc_target;
setappdata(ax, 'viz_last_xc', xc);

XL = [xc - halfW, xc + halfW];
YL = [0, 30];   % Fixed road Y-range [m]

% Set limits without calling axis(ax,'manual') to avoid listener crashes.
set(ax, 'XLim',XL, 'YLim',YL);

%% ── Road Background Bands ────────────────────────────────────────────────
%   Two shaded bands represent the two carriageways of the highway.
%   Bands are pushed to the bottom of the rendering stack.

bands     = [3 13; 15 25];
bandFace  = [0.92 0.92 0.92];
bandAlpha = 0.8;

for b = 1:size(bands, 1)
    xPatch = [XL(1)-150, XL(2)+150, XL(2)+150, XL(1)-150];
    yPatch = [bands(b,1), bands(b,1), bands(b,2), bands(b,2)];
    hBand  = patch(ax, xPatch, yPatch, bandFace, ...
        'FaceAlpha',bandAlpha, 'EdgeColor','none', ...
        'HandleVisibility','off', 'Clipping','on');
    uistack(hBand, 'bottom');
end

%% ── Title and Axis Labels ────────────────────────────────────────────────
ttl = sprintf('$T_\\mathrm{S}$=%.2f s | %d vehicles | %d interferer paths', ...
    t0, height(S), height(R));
if ~isempty(vehIdOpt)
    ttl = sprintf('%s | Victim Vehicle ID: %g', ttl, vehIdOpt);
end
title(ax,  ttl);
xlabel(ax, 'Longitudinal distance $x$ (m)', 'Interpreter','latex');
ylabel(ax, 'Lateral distance $y$ (m)',       'Interpreter','latex');

%% ── Speed Colormap Setup ─────────────────────────────────────────────────
%   Colorbar is hard-clamped to [100 160] km/h to match the highway
%   scenario speed range used in the published case study.
%   cfun maps a speed value to its RGB colour for vehicle rendering.

smin = prctile(S.V_speed, 5);
smax = prctile(S.V_speed, 95);
if ~isfinite(smin) || ~isfinite(smax) || smax <= smin
    smin = min(S.V_speed);
    smax = max(S.V_speed) + eps;
end

cfun = @(v) cmap(max(1, min(256, floor(1 + 255*(v - smin) / max(eps, smax - smin)))), :);

clim(ax, [100 160]);
set(hcb, 'Limits', [100 160]);
hcb.Label.String          = 'Vehicle speed (km/h)';
hcb.Label.Interpreter     = 'latex';
hcb.TickLabelInterpreter  = 'latex';

%% ── Lane Assignment ──────────────────────────────────────────────────────
%   Lane centres are estimated via k-means on the first call and cached
%   in appdata.  This prevents lane assignment from flickering between
%   frames due to k-means stochasticity.

laneEdges = [];
if isappdata(ax, 'viz_lane_edges')
    laneEdges = getappdata(ax, 'viz_lane_edges');
end

if isempty(laneEdges) || numel(laneEdges) ~= NumLanes + 1
    yall = double(S.V_y);
    try
        ctrs      = sort(kmeans(yall, NumLanes, 'Replicates',5, 'MaxIter',200));
        laneEdges = [-inf, (ctrs(1:end-1) + ctrs(2:end))/2, inf];
    catch
        qs        = quantile(yall, linspace(0, 1, NumLanes+1));
        laneEdges = [-inf, qs(2:end-1), inf];
    end
    setappdata(ax, 'viz_lane_edges', laneEdges);
end

laneIdx = discretize(double(S.V_y), laneEdges);

%% ── Vehicle Rendering ────────────────────────────────────────────────────
%   Each vehicle is rendered as either a chamfered car icon (useCars = true)
%   or a raw corner polygon (useCars = false).
%   Heading angle is determined by lane direction (east or west).

for k = 1:height(S)

    xpoly = [S.corner_p1_x(k), S.corner_p2_x(k), S.corner_p3_x(k), S.corner_p4_x(k)];
    ypoly = [S.corner_p1_y(k), S.corner_p2_y(k), S.corner_p3_y(k), S.corner_p4_y(k)];

    % Fall back to centred rectangle if corner coordinates are missing.
    if any(~isfinite(xpoly)) || any(~isfinite(ypoly))
        xc_veh = double(S.V_x(k));
        yc_veh = double(S.V_y(k));
        Lveh   = mean(CarLenClamp);
        Wveh   = mean(CarWidClamp);
    else
        xc_veh = mean(xpoly);
        yc_veh = mean(ypoly);
        Lveh   = max(xpoly) - min(xpoly);
        Wveh   = max(ypoly) - min(ypoly);
        if ~isfinite(Lveh) || Lveh <= 0,  Lveh = mean(CarLenClamp);  end
        if ~isfinite(Wveh) || Wveh <= 0,  Wveh = mean(CarWidClamp);  end
    end

    % Clamp vehicle dimensions to configured bounds.
    Lveh = min(max(Lveh, CarLenClamp(1)), CarLenClamp(2));
    Wveh = min(max(Wveh, CarWidClamp(1)), CarWidClamp(2));

    % Assign heading from lane direction.
    li = laneIdx(k);
    if isnan(li),  li = 1;  end
    theta = 180 * ismember(li, LanesWest);   % 180° = westbound, 0° = eastbound

    if useCars
        draw_car_icon(ax, xc_veh, yc_veh, Lveh, Wveh, theta, cfun(S.V_speed(k)));
    else
        if any(~isfinite(xpoly)) || any(~isfinite(ypoly))
            dx    = Lveh/2;  dy = Wveh/2;
            xpoly = [xc_veh-dx, xc_veh+dx, xc_veh+dx, xc_veh-dx];
            ypoly = [yc_veh-dy, yc_veh-dy, yc_veh+dy, yc_veh+dy];
        end
        patch('Parent',ax, 'XData',xpoly, 'YData',ypoly, ...
            'FaceColor',cfun(S.V_speed(k)), 'FaceAlpha',0.95, ...
            'EdgeColor',[0.2 0.2 0.2], 'LineWidth',1.0, ...
            'HandleVisibility','off');
    end

    % Optional heading debug arrows (first 8 vehicles only).
    if dbgHead && k <= 8
        quiver(ax, xc_veh, yc_veh, cosd(theta), sind(theta), 5, 'k', ...
            'MaxHeadSize',1, 'HandleVisibility','off');
    end

end   % vehicle loop

%% ── Interference Path Downselection ─────────────────────────────────────
%   When the number of paths exceeds maxPaths, retain only the longest
%   paths (by total LOS + reflected path length) to prioritise the most
%   severe interferers visually.

if height(R) > maxPaths
    totLen  = fillnan(R.lengthPathFirst) + fillnan(R.lengthPathSecond);
    [~, ord] = sort(totLen, 'descend');
    R        = R(ord(1:maxPaths), :);
end

%% ── Interference Path Rendering ──────────────────────────────────────────
%   LOS paths: solid red line between tx and rx sensor positions.
%   Reflected paths: two dashed grey segments via the reflection point,
%                    with a yellow marker at the reflection point.

for i = 1:height(R)

    tx = R.interferingSensor_x(i);   ty = R.interferingSensor_y(i);
    rx = R.interferedSensor_x(i);    ry = R.interferedSensor_y(i);
    px = R.reflectPointX(i);         py = R.reflectPointY(i);
    L2 = R.lengthPathSecond(i);

    % Classify as LOS if reflection point or second-leg length is absent.
    isLOS = ~isfinite(px) || ~isfinite(py) || ~isfinite(L2) || (L2 == 0);

    if isLOS
        plot(ax, [tx rx], [ty ry], '-', ...
            'Color',[0.85 0 0], 'LineWidth',2.0, 'HandleVisibility','off');
    else
        % First leg: transmitter → reflection point
        plot(ax, [tx px], [ty py], '--', ...
            'Color',[0.5 0.5 0.5], 'LineWidth',2.0, 'HandleVisibility','off');
        % Second leg: reflection point → victim receiver
        plot(ax, [px rx], [py ry], '--', ...
            'Color',[0.5 0.5 0.5], 'LineWidth',2.0, 'HandleVisibility','off');
        % Reflection point marker
        plot(ax, px, py, 'o', 'MarkerSize',6, ...
            'MarkerFaceColor','y', 'MarkerEdgeColor','k', ...
            'LineWidth',0.8, 'HandleVisibility','off');
    end

end   % path loop

%% ── Focus Vehicle Highlight ──────────────────────────────────────────────
%   A blue circle is drawn around the focus vehicle's bounding box centre.
%   The circle radius is proportional to the vehicle's bounding diagonal.

if ~isempty(vehIdOpt)
    r = S(S.Vid == vehIdOpt, :);
    if ~isempty(r)
        xpoly = [r.corner_p1_x, r.corner_p2_x, r.corner_p3_x, r.corner_p4_x];
        ypoly = [r.corner_p1_y, r.corner_p2_y, r.corner_p3_y, r.corner_p4_y];

        if all(isfinite(xpoly)) && all(isfinite(ypoly))
            xc0   = mean(xpoly);
            yc0   = mean(ypoly);
            dx    = max(xpoly) - min(xpoly);
            dy    = max(ypoly) - min(ypoly);
            rCirc = 0.8 * hypot(dx, dy);   % radius proportional to bounding diagonal

            th = linspace(0, 2*pi, 180);
            hC = plot(ax, xc0 + rCirc*cos(th), yc0 + rCirc*sin(th), ...
                'b-', 'LineWidth',3, 'HandleVisibility','off', 'Clipping','on');
            uistack(hC, 'top');
        end
    end
end

% Reassert axis limits after all rendering (prevents auto-rescaling).
set(ax, 'XLim',XL, 'YLim',YL);

end   % ARIFOMA_draw_scene_snapshot

% =========================================================================
%% ── Local Helper Functions ───────────────────────────────────────────────
% =========================================================================

% -------------------------------------------------------------------------
function intf = normalize_intf_columns(intf)
% NORMALIZE_INTF_COLUMNS
%   Standardise interference table column names to a canonical set.
%
%   WiLabVISim output files occasionally use inconsistent capitalisation
%   or underscoring in column names across versions.  This function
%   performs a case-insensitive remap to the canonical names used
%   throughout ARIFOMA.
%
%   Canonical mappings (lower-case source → canonical target):
%     interferedsensor_x  → interferedSensor_x
%     interferedsensor_y  → interferedSensor_y
%     interferingsensor_x → interferingSensor_x
%     interferingsensor_y → interferingSensor_y
%     lengthpathfirst     → lengthPathFirst
%     lengthpathsecond    → lengthPathSecond
%     reflectpointx       → reflectPointX
%     reflectpointy       → reflectPointY
%
%   INPUT / OUTPUT
%     intf  (table) — interference table; column names are updated in-place.

    canon = { ...
        'interferedsensor_x',   'interferedSensor_x';
        'interferedsensor_y',   'interferedSensor_y';
        'interferingsensor_x',  'interferingSensor_x';
        'interferingsensor_y',  'interferingSensor_y';
        'lengthpathfirst',      'lengthPathFirst';
        'lengthpathsecond',     'lengthPathSecond';
        'reflectpointx',        'reflectPointX';
        'reflectpointy',        'reflectPointY'};

    names = intf.Properties.VariableNames;
    L     = lower(names);

    for k = 1:size(canon, 1)
        idx = find(strcmp(L, canon{k,1}), 1);
        if ~isempty(idx)
            names{idx} = canon{k,2};
            L{idx}     = lower(canon{k,2});
        end
    end

    intf.Properties.VariableNames = names;
end

% -------------------------------------------------------------------------
function v = fillnan(v)
% FILLNAN  Replace non-finite values with zero.
%
%   Used to safely sum LOS and reflected path lengths when either leg may
%   be NaN (e.g. LOS-only paths have no second leg).

    v(~isfinite(v)) = 0;
end

% -------------------------------------------------------------------------
function draw_car_icon(ax, xc, yc, L, W, heading_deg, faceColor)
% DRAW_CAR_ICON  Render a chamfered rectangle car glyph at a given pose.
%
%   Produces a top-down car silhouette as a MATLAB patch object.
%   The body shape is a slightly chamfered rectangle (clipped corners)
%   that reads more naturally as a vehicle than a plain rectangle.
%
%   INPUTS
%     ax           — target axes handle
%     xc, yc       — centre position [m]
%     L            — vehicle length [m]
%     W            — vehicle width [m]
%     heading_deg  — heading angle [degrees]; 0° = east, 180° = west
%     faceColor    — RGB face colour [1×3 double]

    th  = deg2rad(heading_deg);
    R   = [cos(th), -sin(th); sin(th), cos(th)];
    rot = @(X, Y) (R * [X(:)'; Y(:)'])';

    % Chamfered rectangle body: 8 vertices with clipped corners.
    r  = 0.20 * min(L, W);   % chamfer radius proportional to smallest dimension
    xb = [-L/2+r,  L/2-r,  L/2,    L/2,    L/2-r, -L/2+r, -L/2,   -L/2  ];
    yb = [-W/2,   -W/2,   -W/2+r,  W/2-r,  W/2,    W/2,    W/2-r, -W/2+r];

    B = rot(xb, yb) + [xc, yc];

    patch('Parent',ax, 'XData',B(:,1), 'YData',B(:,2), ...
        'FaceColor',faceColor, 'FaceAlpha',0.95, ...
        'EdgeColor',[0.2 0.2 0.2], 'LineWidth',0.9, ...
        'HandleVisibility','off');

end   % draw_car_icon