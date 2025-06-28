close all; clc; clearvars -except ps5000aDeviceObj;

%% === Experiment Parameters ===       
f = 4.8e6;                          % Frequency of the piezo
T = 1/f;                            % Period of the square wave
t_window = 0.5e-6;                    % Set equal to pulse width from function generator
N = t_window/T;                     % Number of square wave periods in each pulse
PreTrig = 0;                        % Number of samples to be gathered from before the trigger
PostTrig = 3000;                    % Number of samples to be gathered after the trigger
M = 2^4;                           % Average over M acquired signals to reduce noise
bandpass = [3e6 7e6];               % Set frequnecy range for the bandpass filter
order = 4;                          % Order of the bandpass filter
%% === Stage setup and GUI ===
global hx hy hz
f1 = figure('Position', [20, 502, 600, 450], 'Menu', 'None', 'Name', 'APT GUI - X Axis');
f2 = figure('Position', [650, 502, 600, 450], 'Menu', 'None', 'Name', 'APT GUI - Y Axis');
f3 = figure('Position', [650, 20, 600, 450], 'Menu', 'None', 'Name', 'APT GUI - Z Axis');
hx = actxcontrol('MGMOTOR.MGMotorCtrl.1', [20 20 600 400], f1);
hy = actxcontrol('MGMOTOR.MGMotorCtrl.1', [20 20 600 400], f2);
hz = actxcontrol('MGMOTOR.MGMotorCtrl.1', [20 20 600 400], f3);
hx.StartCtrl; hy.StartCtrl; hz.StartCtrl;
set(hx,'HWSerialNum',27266153); set(hy,'HWSerialNum',27266066); set(hz,'HWSerialNum',27266157);
hx.registerevent({'MoveComplete' 'MoveCompleteHandler'});
hy.registerevent({'MoveComplete' 'MoveCompleteHandler'});
hz.registerevent({'MoveComplete' 'MoveCompleteHandler'});
hx.SetVelParams(0,1,4.5,2.4); hy.SetVelParams(0,1,4.5,2.4); hz.SetVelParams(0,1,4.5,2.4);

%% === PicoScope Initialization ===
if exist('ps5000aDeviceObj','var') && isvalid(ps5000aDeviceObj)
    warning('Previous PicoScope session detected. Disconnecting...')
    disconnect(ps5000aDeviceObj);
    clear ps5000aDeviceObj
end
PS5000aConfig;
[ps5000aStructs, ps5000aEnuminfo] = ps5000aSetConfig();
ps5000aDeviceObj = icdevice('picotech_ps5000a_generic','');
connect(ps5000aDeviceObj);

%% === Channel Setup === 
% Channel A
chA = ps5000aEnuminfo.enPS5000AChannel.PS5000A_CHANNEL_A;
channelSettings(1).enabled = PicoConstants.TRUE;
channelSettings(1).coupling = ps5000aEnuminfo.enPS5000ACoupling.PS5000A_DC;
channelSettings(1).range = ps5000aEnuminfo.enPS5000ARange.PS5000A_10V;
channelSettings(1).analogueOffset = 0.0;
%channelARangeMv = PicoConstants.SCOPE_INPUT_RANGES(channelSettings(1).range + 1);
% Channel B
chB = ps5000aEnuminfo.enPS5000AChannel.PS5000A_CHANNEL_B;
channelSettings(2).enabled = PicoConstants.TRUE;
channelSettings(2).coupling = ps5000aEnuminfo.enPS5000ACoupling.PS5000A_DC;
channelSettings(2).range = ps5000aEnuminfo.enPS5000ARange.PS5000A_500MV;
channelSettings(2).analogueOffset = 0.0;
%channelBRangeMv = PicoConstants.SCOPE_INPUT_RANGES(channelSettings(2).range + 1);

% Keep the status values returned from the driver.
numChannels = get(ps5000aDeviceObj, 'channelCount');
status.setChannelStatus = zeros(numChannels, 1);
[status.currentPowerSource] = invoke(ps5000aDeviceObj, 'ps5000aCurrentPowerSource');
% Check if the power supply is connected - channels C and D will not be
% enabled on a 4-channel oscilloscope if it is only USB powered.
if (status.currentPowerSource == PicoStatus.PICO_POWER_SUPPLY_NOT_CONNECTED)
    numChannels = PicoConstants.DUAL_SCOPE;
end
for ch = 1:numChannels
    status.setChannelStatus(ch) = invoke(ps5000aDeviceObj, 'ps5000aSetChannel', ...
    (ch - 1), channelSettings(ch).enabled, ...
    channelSettings(ch).coupling, channelSettings(ch).range, ...
    channelSettings(ch).analogueOffset);
end

%% === Set Resolution and Timebase
[status.setResolution, scope.resolution] = invoke(ps5000aDeviceObj, 'ps5000aSetDeviceResolution', 15);
maxADCCount = get(ps5000aDeviceObj, 'maxADCValue');

%Timebase determines available samples in the selected segment
scope.timebaseIndex = 3;
status.getTimebase = PicoStatus.PICO_INVALID_TIMEBASE;
while (status.getTimebase == PicoStatus.PICO_INVALID_TIMEBASE)
[status.getTimebase, scope.timeIntervalNanoSeconds, scope.maxSamples] = invoke(ps5000aDeviceObj, 'ps5000aGetTimebase', scope.timebaseIndex, 0);
if (status.getTimebase == PicoStatus.PICO_OK)
break;
else
scope.timebaseIndex = scope.timebaseIndex + 1;
end
end
set(ps5000aDeviceObj, 'timebase', scope.timebaseIndex);

fs = 125e6/(scope.timebaseIndex-2);


%% === Set center of measurement plane === 
%[x_0, y_0, z_0] = move_stage_gui(hx, hy, hz);

%Option with live plotting
%[x_0, y_0, z_0, filename] = move_stage_gui_live(hx, hy, hz, ps5000aDeviceObj, scope.timeIntervalNanoSeconds, bandpass, N, f);

%Option with live plotting and range
[x_0, y_0, z_0, filename, sensitivity] = move_stage_gui_live_wRange(hx, hy, hz, ps5000aDeviceObj, scope.timeIntervalNanoSeconds, bandpass, N, f);

%% === Set scan plane and step size ===
%[ax1, ax2, fixed_axis, ax1_start, ax2_start, fixed_val, plane, grid_size, grid_step, axis1Vals, axis2Vals, fixedVal] = select_scan_plane_with_motors(hx, hy, hz);
[ax1, ax2, fixed_axis, ax1_start, ax2_start, fixed_val, plane, grid_size, grid_step, axis1Vals, axis2Vals, fixedVal] = select_scan_plane_with_motors_rectangle(hx, hy, hz);

%% === Trigger Setup === 
%A trigger on channel A will be used to indicate when data should be gathered
triggerGroupObj = get(ps5000aDeviceObj, 'Trigger');
triggerGroupObj = triggerGroupObj(1);
set(triggerGroupObj, 'autoTriggerMs', 10000); %Automatically gather after x ms without trigger

Channel = 0; %Channel A (ps5000aEnuminfo.enPS5000AChannel.PS5000A_CHANNEL_A)
Threshold = 1000; %mV
Direction = 2; %Rising
[status.setSimpleTrigger] = invoke(triggerGroupObj, 'setSimpleTrigger', Channel, Threshold, Direction);

%% === Rapid Block Mode Configuration
[status.memorySegments, nMaxSamples] = invoke(ps5000aDeviceObj, 'ps5000aMemorySegments', M);
set(ps5000aDeviceObj, 'numPreTriggerSamples', PreTrig);
set(ps5000aDeviceObj, 'numPostTriggerSamples', PostTrig);

rapidBlockGroupObj = get(ps5000aDeviceObj, 'Rapidblock'); rapidBlockGroupObj = rapidBlockGroupObj(1);
invoke(rapidBlockGroupObj, 'ps5000aSetNoOfCaptures', M);
blockObj = get(ps5000aDeviceObj, 'Block'); blockObj = blockObj(1);


%% === Initialize result structure
Nx = length(axis1Vals);
Ny = length(axis2Vals);
results(Nx, Ny) = struct('x', [], 'y', [], 'z', [], 'vpp', [], 'raw', []);

%% === Raster scan loop === 
totalPoints = Nx*Ny;

for ix = 1:Nx
    for iy = 1:Ny
        moveStageTo(ax1, axis1Vals(ix));
        if iy == 1
            pause(1)
        end
        moveStageTo(ax2, axis2Vals(iy));
        moveStageTo(fixed_axis, fixedVal);
        pause(0.2)

        invoke(blockObj, 'runBlock', 0);
        isReady = false;
        while ~isReady
            [~, isReady] = invoke(blockObj, 'ps5000aIsReady');
            pause(0.01);
        end
        [~, ~, chA_all, chB_all] = invoke(rapidBlockGroupObj, 'getRapidBlockData', M, 1, 0);
        avg_signal = mean(chB_all, 2)/2;
        t = (0: PreTrig + PostTrig - 1) * double(scope.timeIntervalNanoSeconds) / 1e6;

        burstEnd = floor(N/f*fs);
        signalSegment = avg_signal;
        signalSegment(1:burstEnd) = [];
        [b,a] = butter(order, bandpass/(fs/2), 'bandpass');
        filtered = filtfilt(b, a, signalSegment);

        vpp = max(filtered) - min(filtered);

        results(ix,iy).x = axis1Vals(ix);
        results(ix,iy).y = axis2Vals(iy);
        results(ix,iy).z = fixedVal;
        results(ix,iy).vpp = vpp;
        results(ix,iy).raw = avg_signal;
        results(ix,iy).rawpres = avg_signal/sensitivity;
        results(ix,iy).pppres = vpp/sensitivity;
        
        figure(99); clf;
        subplot(2,1,1); plot(t, avg_signal/sensitivity); title('Raw Signal'); xlabel('Time (ms)'); ylabel('pressure (MPa)');
        subplot(2,1,2); plot(t(burstEnd+1:end), filtered/sensitivity); title(sprintf('Filtered | Pressure =%.2f MPa', vpp/sensitivity)); ylabel('Pressure (MPa)');
        drawnow

        fprintf('(%d, %d) @ [%.2f, %.2f, %.2f], Vpp = %.2f mV\n', ix, iy, axis1Vals(ix), axis2Vals(iy), fixedVal, vpp);
        currentPoint = (ix-1) * Ny + iy;
        fprintf('%3d/%3d', currentPoint, totalPoints)
    end
    %save(filename);
end
save(filename);
%% === Return to start coordinates
moveStageTo(ax1, ax1_start);
moveStageTo(ax2, ax2_start);
moveStageTo(fixed_axis, fixed_val);

%% === Colourmap of peak-to-peak pressure === 
vpp_map = arrayfun(@(results) results.pppres, results');
figure;
imagesc(axis1Vals, axis2Vals, vpp_map);
axis xy;
colorbar;
xlabel(sprintf('%s (mm)', plane(1)));
ylabel(sprintf('%s (mm)', plane(2)));
title('Pressure Distribution (MPa)')


%% === Disconnect from Picoscope === 
disconnect(ps5000aDeviceObj);
clear ps5000aDeviceObj

%% === Helper Function Movement === 
function moveStageTo(motor, pos)
    motor.SetAbsMovePos(0, pos);
    motor.MoveAbsolute(0, false);
    status = motor.GetStatusBits_Bits(0);
    while IsMoving(status)
        status = motor.GetStatusBits_Bits(0);
        pause(0.05);
    end
end

