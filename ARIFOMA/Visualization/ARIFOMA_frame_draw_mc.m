function ARIFOMA_frame_draw_mc(axTop, frame, opt)
% ARIFOMA_FRAME_DRAW_MC
%   Render a pre-extracted scene frame into a target axes handle.
%
% PURPOSE
% -------
%   This function is a pure renderer — it accepts a frame struct that has
%   already been extracted by ARIFOMA_frame_extract_mc and delegates all
%   drawing to ARIFOMA_draw_scene_snapshot.
%
%   The separation between extraction (ARIFOMA_frame_extract_mc) and
%   rendering (this function) allows the same frame to be redrawn multiple
%   times (e.g. for different axes or export resolutions) without
%   re-querying the scenario tables.
%
% INPUTS
%   axTop  — target MATLAB axes handle to draw into
%   frame  — scene frame struct produced by ARIFOMA_frame_extract_mc;
%            required fields:
%              frame.Vin    — vehicle location table sliced at frame.t0
%              frame.R_draw — interference path table sliced at frame.t0
%              frame.t0     — actual snapshot timestamp [s]
%              frame.vehId  — focus vehicle ID (victim)
%              frame.XL     — optional x-axis limits [xMin xMax];
%                             empty [] disables the override
%   opt    — configuration struct with rendering knobs; required fields:
%              opt.Tol      — time tolerance for snapshot slicing [s]
%              opt.MaxPaths — maximum number of interference paths to draw
%              opt.StretchY — vertical display stretch factor [0.1, 1]
%
% OUTPUTS
%   None — all output is rendered directly into axTop.
%
% CONTRACT
% --------
%   This function must not modify frame, opt, or any simulation state.
%   It is safe to call repeatedly inside a frame loop.
%
% SEE ALSO
%   ARIFOMA_frame_extract_mc, ARIFOMA_draw_scene_snapshot, main_sim

% =========================================================================

%% ── Axes Preparation ─────────────────────────────────────────────────────
%   Set axTop as the current axes and clear all existing graphics objects
%   to ensure a clean slate before rendering the new frame.
%   The colorbar cached by ARIFOMA_draw_scene_snapshot survives this
%   deletion because it is stored in appdata, not as a child of axTop.

axes(axTop);
delete(allchild(axTop));

%% ── Scene Rendering ──────────────────────────────────────────────────────
%   Delegate all drawing to ARIFOMA_draw_scene_snapshot.
%   The focus vehicle (frame.vehId) is passed so that its interference
%   paths are highlighted and the camera follows it across frames.

ARIFOMA_draw_scene_snapshot( ...
    frame.Vin, ...
    frame.R_draw, ...
    frame.t0, ...
    opt.Tol, ...
    opt.MaxPaths, ...
    opt.StretchY, ...
    axTop, ...
    'VehicleId', frame.vehId, ...
    'WinX',      100);

end   % ARIFOMA_frame_draw_mc