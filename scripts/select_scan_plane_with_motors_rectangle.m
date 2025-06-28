function [ax1, ax2, fixed_axis, ax1_loc, ax2_loc, fixed_loc, plane, ...
grid_size, grid_step, axis1Vals, axis2Vals, fixedVal] = ...
select_scan_plane_with_motors_rectangle(hx, hy, hz)

%% OUTPUT INITIALIZATION

safeRange.x = [0, 21];
safeRange.y = [0, 26];
safeRange.z = [20, 50];

ax1 = []; ax2 = []; fixed_axis = [];
ax1_loc = []; ax2_loc = []; fixed_loc = [];
plane = ''; grid_size = []; grid_step = [];
axis1Vals = []; axis2Vals = []; fixedVal = [];

% Create GUI
fig = uifigure('Name','Scan Plane Selection','Position',[100 100 900 450]);

% Plane selection
uilabel(fig,'Text','Scan plane:', 'Position',[20 400 100 22]);
planeDrop = uidropdown(fig, ...
'Items',{'XY','XZ','YZ'}, ...
'Position',[120 400 100 22], ...
'ValueChangedFcn', @(~,~)updatePreview());

% Grid size and step input
uilabel(fig, 'Text', 'Size ax1 (mm):', 'Position', [20 360 100 20]);
size1Field = uieditfield(fig, 'numeric', 'Position', [130 360 80 22], ...
    'Value', 2.0, 'ValueChangedFcn', @(~,~)updatePreview());
uilabel(fig, 'Text', 'Size ax2 (mm):', 'Position', [20 330 100 20]);
size2Field = uieditfield(fig, 'numeric', 'Position', [130 330 80 22], ...
    'Value', 2.0, 'ValueChangedFcn', @(~,~)updatePreview());
uilabel(fig, 'Text', 'Grid step (mm):', 'Position', [20 300 100 20]);
stepField = uieditfield(fig, 'numeric', 'Position', [130 300 80 22], ...
'Value', 0.2, 'ValueChangedFcn', @(~,~)updatePreview());

% Label to show estimated number of points
pointLabel = uilabel(fig, 'Text', '', 'Position', [20 280 200 22]);

% Confirm button
btnConfirm = uibutton(fig, 'Text', 'Confirm and continue', ...
    'Position', [30 220 180 30], ...
    'ButtonPushedFcn', @(~,~)finalizeSelection());

% Warning label
warningLabel = uilabel(fig, 'Text', '', ...
    'Position', [30 190 250 22], 'FontColor', [0.8 0 0]);

% 3D preview axis
ax = uiaxes(fig, 'Position', [300 80 550 320]);
xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z'); grid(ax,'on');
view(ax, 3);
axis(ax, [0 30 0 30 0 50]);

% Draw safe range cube as dashed outline
hold(ax, 'on');
plot3(ax, ...
    [safeRange.x(1), safeRange.x(2), safeRange.x(2), safeRange.x(1), safeRange.x(1), ...
     safeRange.x(1), safeRange.x(2), safeRange.x(2), safeRange.x(1), safeRange.x(1), ...
     safeRange.x(1), safeRange.x(1)], ...
    [safeRange.y(1), safeRange.y(1), safeRange.y(2), safeRange.y(2), safeRange.y(1), ...
     safeRange.y(1), safeRange.y(1), safeRange.y(2), safeRange.y(2), safeRange.y(2), ...
     safeRange.y(2), safeRange.y(1)], ...
    [safeRange.z(1), safeRange.z(1), safeRange.z(1), safeRange.z(1), safeRange.z(1), ...
     safeRange.z(2), safeRange.z(2), safeRange.z(2), safeRange.z(2), safeRange.z(1), ...
     safeRange.z(2), safeRange.z(2)], ...
    'Color', [0.5 0.5 0.5], 'LineStyle', '--');
hold(ax, 'off');

% Patch object for scan plane
h_scan = patch(ax, 'Faces', [], 'Vertices', [], ...
'FaceColor', [0 0.6 0], 'FaceAlpha', 0.4, 'EdgeColor', 'k');

%% Launch GUI and wait for user input
updatePreview();
uiwait(fig);

%% NESTED FUNCTIONS

function updatePreview()
    % Get user input
    p = planeDrop.Value;
    sz1 = size1Field.Value;
    sz2 = size2Field.Value;
    st = stepField.Value;
    if st <= 0 || sz1 <= 0 || sz2 <= 0
        pointLabel.Text = 'Invalid input.';
        return
    end

    % Determine axes
    switch p
        case 'XY', ax1 = hx; ax2 = hy; fixed_axis = hz;
        case 'XZ', ax1 = hx; ax2 = hz; fixed_axis = hy;
        case 'YZ', ax1 = hy; ax2 = hz; fixed_axis = hx;
    end

    % Get current motor positions
    ax1_loc = ax1.GetPosition_Position(0);
    ax2_loc = ax2.GetPosition_Position(0);
    fixed_loc = fixed_axis.GetPosition_Position(0);
    fixedVal = fixed_loc;

    % Compute scan grid
    range1 = sz1 / 2;
    range2 = sz2 / 2;
    axis1Vals = ax1_loc - range1 : st : ax1_loc + range1;
    axis2Vals = ax2_loc - range2 : st : ax2_loc + range2;
    N1 = numel(axis1Vals);
    N2 = numel(axis2Vals);
    pointLabel.Text = sprintf('Estimated: %d × %d points (%d total)', N1, N2, N1*N2);

    % Create preview plane (4 corner points)
    [X, Y] = meshgrid([min(axis1Vals), max(axis1Vals)], ...
                      [min(axis2Vals), max(axis2Vals)]);
    switch p
        case 'XY'
            verts = [X(:), Y(:), fixed_loc * ones(4,1)];
        case 'XZ'
            verts = [X(:), fixed_loc * ones(4,1), Y(:)];
        case 'YZ'
            verts = [fixed_loc * ones(4,1), X(:), Y(:)];
    end
    h_scan.Vertices = verts;
    h_scan.Faces = [1 2 4 3];

    % Check whether all corners are within safe range
    allX = verts(:,1); allY = verts(:,2); allZ = verts(:,3);
    inX = all(allX >= safeRange.x(1) & allX <= safeRange.x(2));
    inY = all(allY >= safeRange.y(1) & allY <= safeRange.y(2));
    inZ = all(allZ >= safeRange.z(1) & allZ <= safeRange.z(2));

    if inX && inY && inZ
        h_scan.FaceColor = [0 0.6 0]; % safe = green
        btnConfirm.Enable = true;
        warningLabel.Text = '';
    else
        h_scan.FaceColor = [0.8 0 0]; % out of bounds = red
        btnConfirm.Enable = false;
        warningLabel.Text = '⚠️ Scan plane exceeds safe range.';
    end
end

function finalizeSelection()
    % Prevent exit if scan range is invalid (failsafe)
    if ~btnConfirm.Enable
        uialert(fig, 'Scan plane is outside the safe range. Please adjust before continuing.', 'Invalid Selection');
        return;
    end
    % Store output values
    plane = planeDrop.Value;
    grid_size = [size1Field.Value, size2Field.Value];
    grid_step = stepField.Value;
    uiresume(fig);
    delete(fig);
end
end