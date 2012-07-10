


% emulate localization performance with ideal time series data
clear
close all
clc

% sim parameters
fs = 236660.3;                   % sampling rate of generated waveform [Hz]
srcLoc = [randn(1,2)/5 1.5];      % (x,y,z) coordinates
N0 = 10^(-6/20);          % add uncorrelated noise

% create time series of emitted pulse
f0 = 105e3;             % initial frequency [Hz]
f1 = 15e3;              % final frequency [Hz]
T = 0.001;              % pulse length [sec]

t = (0:1/fs:T)';
IF = f0*ones(size(t)) + (f1-f0)/(T).*t;
IA = raisedcos(length(IF));
ts.time = t;
ts.data = real(gen_ifpulse(fs,IF,0,IA));
ts.fs = fs;

% calculate array positions
array = az_positions([6 6],[.3048 .254]);
micLoc = [array.xPos' array.yPos' zeros(length(array.xPos),1)];

%% apply different array perturbations
%micLoc = [micLoc ; micLoc + ones(9,1)*[0 0 1]];
%micLoc = [micLoc ; [0 0 .001]];
%micLoc(5,3) = 1e-3 * randn(1,1);
micLoc(:,3) = micLoc(:,3) + 1e-3 * randn(length(array.xPos),1);

% rotate matrix about x axis
if 1
    theta = -0*pi/180;
    Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
    %srcLoc = (Rx * srcLoc')';   % rotate source
    for k = 1:size(micLoc,1)
        micLoc(k,:) = (Rx * micLoc(k,:)')';     % rotate mics
    end
end

% simulate appropriate delay to each channel
sim = array_sim(ts, micLoc, srcLoc);

L = size(sim.data,1);
M = size(sim.data,2);

% add uncorrelated noise
sim.data = sim.data + N0 * randn(L,M);

% localize using TDOA
coords = TDOA_frame(sim.data', micLoc, sim.fs);

% show results
fprintf('Calculated xyz position: (%g,  %g,  %g)\n',coords(1:3));
fprintf('Actual xyz position:     (%g,  %g,  %g)\n',srcLoc);
fprintf('XYZ error:               (%g,  %g,  %g)\n',coords(1:3) - srcLoc);

% plot spatial representation with source locations (Note: Y and Z are reversed to simplify 3d rotation)
figure
plot3(micLoc(:,1),micLoc(:,3),micLoc(:,2),'.')      % sensor locations
grid on
hold on
set(gca,'YDir','reverse')
xlabel('x [m]')
ylabel('z [m]')
zlabel('y [m]')
plot3(srcLoc(:,1),srcLoc(:,3),srcLoc(:,2),'r.')     % actual location
plot3(coords(1),coords(3),coords(2),'ro')           % estimated location
axis equal

% plot time series of each received signal
figure
plot(sim.time, sim.data + 2*ones(L,1)*(0:M-1))
