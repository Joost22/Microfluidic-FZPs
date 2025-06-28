function [x_start, y_start, z_start, filename, sensitivity] = move_stage_gui_live_wRange(hx, hy, hz, ps5000aDeviceObj, timeIntervalNanoSeconds, bandpass, N, f)
% MOVE_STAGE_GUI_LIVE_WRANGE - GUI for manual stage movement with optional live plot, Vpp feedback, and safety range enforcement

x_start = []; y_start = []; z_start = []; sensitivity = [];
signal = []; t = []; fs = []; M2 = []; vpp = [];
filenameIsSet = false;
xMode = "time"; yMode = "voltage"; sensitivity = 179;

% Define safe range limits (in mm)
safeRange.x = [0, 21];
safeRange.y = [0, 26];
safeRange.z = [0, 50];

% GUI window
fig = uifigure('Name', 'Move Stage (with live Vpp)', 'Position', [100, 100, 950, 500]);

% Filename controls
uilabel(fig, 'Position', [30 20 100 22], 'Text', 'File name');
filenameField = uieditfield(fig, 'text', 'Position', [120 20 220 25], 'Value', '');
setFilenameBtn = uibutton(fig, 'Text', 'Set', 'Position', [350 20 60 25], 'ButtonPushedFcn', @(~,~) confirmFileName());

% Step size selection
uilabel(fig, 'Position', [20 450 100 20], 'Text', 'Step size (mm):');
ddStep = uidropdown(fig, 'Items', {'0.01', '0.1', '1'}, 'Position', [120 450 80 22], 'Value', '0.1');

% Motor controls - X
uilabel(fig, 'Position', [20 390 50 20], 'Text', 'X:');
xMinus = uibutton(fig, 'Text', '-', 'Position', [70 390 40 30]);
xPlus = uibutton(fig, 'Text', '+', 'Position', [210 390 40 30]);
xField = uieditfield(fig, 'numeric', 'Position', [120 390 80 30]);
xLED = uilamp(fig, 'Position', [330 395 20 20], 'Color', [0.5 0.5 0.5]);

% Motor controls - Y
uilabel(fig, 'Position', [20 330 50 20], 'Text', 'Y:');
yMinus = uibutton(fig, 'Text', '-', 'Position', [70 330 40 30]);
yPlus = uibutton(fig, 'Text', '+', 'Position', [210 330 40 30]);
yField = uieditfield(fig, 'numeric', 'Position', [120 330 80 30]);
yLED = uilamp(fig, 'Position', [330 335 20 20], 'Color', [0.5 0.5 0.5]);

% Motor controls - Z
uilabel(fig, 'Position', [20 270 50 20], 'Text', 'Z:');
zMinus = uibutton(fig, 'Text', '-', 'Position', [70 270 40 30]);
zPlus = uibutton(fig, 'Text', '+', 'Position', [210 270 40 30]);
zField = uieditfield(fig, 'numeric', 'Position', [120 270 80 30]);
zLED = uilamp(fig, 'Position', [330 275 20 20], 'Color', [0.5 0.5 0.5]);

% Set origin location
setBtn = uibutton(fig, 'Text', 'Set focal spot location', 'Position', [100 200 200 40]);

% Measurement and plotting options
cbPlot = uicheckbox(fig, 'Text', 'Live plot after move', 'Position', [500 450 200 20], 'Value', true);
cbVpp = uicheckbox(fig, 'Text', 'Auto-update Vpp', 'Position', [700 450 200 20], 'Value', true);
uilabel(fig, 'Position', [500 420 100 20], 'Text', 'Averaging:');
editM = uieditfield(fig, 'numeric', 'Position', [600 420 80 22], 'Value', 64);

% Signal plot area
ax = uiaxes(fig, 'Position', [450 80 450 320]);
xlabel(ax, 'Distance (mm)');
ylabel(ax, 'Voltage (mV)');
title(ax, 'Hydrophone signal');
ylim(ax, [-50 50]);

% Axis range and scaling options
uilabel(fig, 'Position', [460 55 80 20], 'Text', 'X-range:');
xRangeField = uieditfield(fig, 'numeric', 'Position', [510 55 50 22], 'Value', 0.04);
xUnitBtn = uibutton(fig, 'Text', 'Time (ms)', 'Position', [570 55 90 22], 'ButtonPushedFcn', @(btn,~) toggleXUnit(btn));
uilabel(fig, 'Position', [730 55 80 20], 'Text', 'Y-range:');
yRangeField = uieditfield(fig, 'numeric', 'Position', [780 55 50 22], 'Value', 100);
yUnitBtn = uibutton(fig, 'Text', 'Voltage (mV)', 'Position', [840 55 90 22], 'ButtonPushedFcn', @(btn,~) toggleYUnit(btn));
uilabel(fig, 'Position', [730 30 80 20], 'Text', 'Sensitivity (mV/MPa):');
sensField = uieditfield(fig, 'numeric', 'Position', [830 30 60 22], 'Value', 179);

% Save signal
uilabel(fig, 'Position', [460 20 100 22], 'Text', 'Tag:');
saveTagField = uieditfield(fig, 'text', 'Position', [510 20 80 25]);
btnSaveSignal = uibutton(fig, 'Text', 'Save signal', 'Position', [600 20 100 25], 'ButtonPushedFcn', @(~,~) saveSignal());

% Button callbacks
xMinus.ButtonPushedFcn = @(~,~) moveAndMeasure(hx, -1, xLED);
xPlus.ButtonPushedFcn = @(~,~) moveAndMeasure(hx, +1, xLED);
yMinus.ButtonPushedFcn = @(~,~) moveAndMeasure(hy, -1, yLED);
yPlus.ButtonPushedFcn = @(~,~) moveAndMeasure(hy, +1, yLED);
zMinus.ButtonPushedFcn = @(~,~) moveAndMeasure(hz, -1, zLED);
zPlus.ButtonPushedFcn = @(~,~) moveAndMeasure(hz, +1, zLED);
xField.ValueChangedFcn = @(src,~) moveTo(hx, src, xLED);
yField.ValueChangedFcn = @(src,~) moveTo(hy, src, yLED);
zField.ValueChangedFcn = @(src,~) moveTo(hz, src, zLED);
setBtn.ButtonPushedFcn = @setStartPosition;

% Apply range updates
xRangeField.ValueChangedFcn = @(src,~) updateAxes();
yRangeField.ValueChangedFcn = @(src,~) updateAxes();
sensField.ValueChangedFcn = @(src,~) updateAxes();

updateFields();
uiwait(fig);

%% ==== Nested functions ====

function moveAndMeasure(motor, direction, lamp)
    step = direction * str2double(ddStep.Value);
    x = hx.GetPosition_Position(0);
    y = hy.GetPosition_Position(0);
    z = hz.GetPosition_Position(0);
    if isequal(motor, hx), x = x + step; end
    if isequal(motor, hy), y = y + step; end
    if isequal(motor, hz), z = z + step; end
    if isWithinSafeRange(x, y, z)
        motor.SetAbsMovePos(0, motor.GetPosition_Position(0) + step);
        motor.MoveAbsolute(0, false);
        waitForMotor(motor, lamp);
        updateFields();
        if cbPlot.Value, doLiveMeasurement(); end
    else
        uialert(fig, '⚠️ Movement blocked: out of safe range.', 'Out of bounds');
    end
end

function moveTo(motor, field, lamp)
    target = round(field.Value, 4);
    x = hx.GetPosition_Position(0);
    y = hy.GetPosition_Position(0);
    z = hz.GetPosition_Position(0);
    if isequal(motor, hx), x = target; end
    if isequal(motor, hy), y = target; end
    if isequal(motor, hz), z = target; end
    if isWithinSafeRange(x, y, z)
        motor.SetAbsMovePos(0, target);
        motor.MoveAbsolute(0, false);
        waitForMotor(motor, lamp);
        updateFields();
        if cbPlot.Value, doLiveMeasurement(); end
    else
        uialert(fig, '⚠️ Movement blocked: out of safe range.', 'Out of bounds');
    end
end

function inBounds = isWithinSafeRange(x, y, z)
    inBounds = ...
        x >= safeRange.x(1) && x <= safeRange.x(2) && ...
        y >= safeRange.y(1) && y <= safeRange.y(2) && ...
        z >= safeRange.z(1) && z <= safeRange.z(2);
end


function waitForMotor(motor, lamp)
    lamp.Color = 'green';
    while IsMoving(motor.GetStatusBits_Bits(0)), pause(0.05); end
    lamp.Color = [0.5 0.5 0.5];
end

function updateFields()
    xField.Value = round(hx.GetPosition_Position(0), 4);
    yField.Value = round(hy.GetPosition_Position(0), 4);
    zField.Value = round(hz.GetPosition_Position(0), 4);
end

function toggleXUnit(btn)
    if strcmp(btn.Text, 'Time (ms)')
        btn.Text = 'Distance (mm)';
    else
        btn.Text = 'Time (ms)';
    end
    updateAxes(); 
end

function toggleYUnit(btn)
    if strcmp(btn.Text, 'Voltage (mV)')
        btn.Text = 'Pressure (MPa)';
    else
        btn.Text = 'Voltage (mV)';
    end
    updateAxes();
end

function setStartPosition(~,~)
    prefix = strtrim(filenameField.Value);
    if isempty(prefix)
        uialert(fig, 'Please enter a prefix before continuing', 'Missing prefix');
        return;
    end
    if ~filenameIsSet
        uialert(fig, 'You must click "Set" to confirm prefix first', 'Prefix not set');
        return; 
    end
    x_start = hx.GetPosition_Position(0);
    y_start = hy.GetPosition_Position(0);
    z_start = hz.GetPosition_Position(0);
    filename = [prefix 'PlaneScan']; 
    uiresume(fig); delete(fig);
end

function doLiveMeasurement()
    M2 = editM.Value;
    PreTrig2 = 0; PostTrig2 = 10000;
    fs = 125e6 / (3 - 2); % Replace with actual timebase logic if needed
    triggerGroupObj = get(ps5000aDeviceObj, 'Trigger'); triggerGroupObj = triggerGroupObj(1);
    set(triggerGroupObj, 'autoTriggerMs', 10000);
    Channel = 0; Threshold = 1000; Direction = 2;
    invoke(triggerGroupObj, 'setSimpleTrigger', Channel, Threshold, Direction);
    [~, ~] = invoke(ps5000aDeviceObj, 'ps5000aMemorySegments', M2);
    set(ps5000aDeviceObj, 'numPreTriggerSamples', PreTrig2);
    set(ps5000aDeviceObj, 'numPostTriggerSamples', PostTrig2);
    rapidBlockGroupObj = get(ps5000aDeviceObj, 'Rapidblock'); rapidBlockGroupObj = rapidBlockGroupObj(1);
    invoke(rapidBlockGroupObj, 'ps5000aSetNoOfCaptures', M2);
    blockObj = get(ps5000aDeviceObj, 'Block'); blockObj = blockObj(1);
    invoke(blockObj, 'runBlock', 0);
    isReady = false;
    while ~isReady
        [~, isReady] = invoke(blockObj, 'ps5000aIsReady');
        pause(0.01);
    end
    [~, ~, ~, chB_all2] = invoke(rapidBlockGroupObj, 'getRapidBlockData', M2, 1, 0);
    signal = mean(chB_all2, 2)/2;
    t = (0:PreTrig2+PostTrig2-1) * double(timeIntervalNanoSeconds) / 1e6;
    if cbVpp.Value
        burstEnd = floor(N/f*fs);
        seg = signal; seg(1:burstEnd) = [];
        [b,a] = butter(4, bandpass/(fs/2), 'bandpass');
        filtered = filtfilt(b, a, seg);
        vpp = max(filtered) - min(filtered);
    end
    % Plot
    x_unit = xUnitBtn.Text;
    y_unit = yUnitBtn.Text;
    x_range = xRangeField.Value;
    y_range = yRangeField.Value;
    sens_mV_per_MPa = sensField.Value;
    sensitivity = sensField.Value;
    if strcmp(x_unit, 'Time (ms)')
        x_data = t;
        xlabel(ax, 'Time (ms)');
    else
        x_data = t * (1481);
        xlabel(ax, 'Distance (mm)');
    end
    if strcmp(y_unit, 'Voltage (mV)')
        y_data = signal;
        ylabel(ax, 'Voltage (mV)');
        ylim(ax, [-y_range/2 y_range/2]);
    else
        y_data = signal / sens_mV_per_MPa;
        ylabel(ax, 'Pressure (MPa)');
        ylim(ax, [-y_range/2 y_range/2]);
    end
    xlim(ax, [0 x_range]);
    cla(ax); plot(ax, x_data, y_data); grid(ax, 'on');
    if cbVpp.Value
        title(ax, sprintf('Live signal: Vpp = %.2f mV', vpp));
    else
        title(ax, 'Live signal');
    end
    drawnow;
end

function updateAxes()
    if isempty(signal) || isempty(t)
        return;
    end

    x_unit = xUnitBtn.Text;
    y_unit = yUnitBtn.Text;
    x_range = xRangeField.Value;
    y_range = yRangeField.Value;
    sens_mV_per_MPa = sensField.Value;
    sensitivity = sensField.Value;

    %X-axis
    if strcmp(x_unit, 'Time (ms)')
        x_data = t;
        xlabel(ax, 'Time (ms)');
    else 
        x_data = t*(1481); 
        xlabel(ax, 'Distance (mm)');
    end

    %Y-axis
    if strcmp(y_unit, 'Voltage (mV)')
        y_data = signal;
        ylabel(ax, 'Voltage (mV)');
    else
        y_data = signal / sens_mV_per_MPa;
        ylabel(ax, 'Pressure (MPa)');
    end

    % Apply limits
    xlim(ax, [0 x_range])
    ylim(ax, [-y_range/2 y_range/2])

    % Update plot
    cla(ax); plot(ax, x_data, y_data); grid(ax, 'on');
    if cbVpp.Value && exist('vpp', 'var')
        title(ax, sprintf('Live signal: Vpp = %.2f', vpp));
    else
        title(ax, 'Live signal');
    end
    drawnow;
end

function saveSignal()
    if ~filenameIsSet
        uialert(fig, 'Please set a valid filename prefix first using the "set" button', 'Prefix is not set');
        return; 
    end
    tag = strtrim(saveTagField.Value);
    prefix = strtrim(filenameField.Value);
    if isempty(tag)
        uialert(fig, 'Please enter a tag befor saving', 'Missing tag');
        return;
    end

    saveName = sprintf('%s_%s.mat', prefix, tag);

    if isfile(saveName)
        uialert(fig,sprintf('File name already exists, choose a different tag'));
        return;
    end

    traceData.signal = signal;
    traceData.vpp = vpp;
    traceData.x = hx.GetPosition_Position(0);
    traceData.y = hy.GetPosition_Position(0);
    traceData.z = hz.GetPosition_Position(0);
    traceData.Sens = sensField.Value;
    traceData.M = editM.Value;
    try 
        save(saveName, 'traceData');
        uialert(fig, sprintf('Signal saved to:\n%s', saveName), 'Success');
    catch ME   
        uialert(fig, sprintf('Error saving:\n%s', ME.message), 'Save Failed');
    end
end

function confirmFileName()
    prefix = strtrim(filenameField.Value);
    if isempty(prefix)
        uialert(fig, 'Please enter a prefix before setting.', 'Missing input');
        return;
    end
    
    mainName = [prefix '_PlaneScan.mat'];
    if isfile(mainName)
        uialert(fig, sprint('A scan already exists for this prefix, choose a different one'), 'Prefix already taken');
        return;
    end

    if ~isempty(dir([prefix '_*.mat']))
        uialert(fig, sprintf('Prefix "%s_" already exists. Choose a unique name.', prefix), 'Conflict');
        return;
    end
    filenameIsSet = true;
    uialert(fig, sprintf('Prefix "%s" is now set', prefix), 'Prefix confirmed');
end


end

