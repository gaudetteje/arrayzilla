function sim = array_sim(ts, sLoc, source)
% ARRAY_SIM  simulates time series signals as received by an array
%
% sim = array_sim(ts, sLoc, source) takes the time series, ts, sensor
%      positions, sLoc, and source location, source, and returns a
%      multi-channel time series, sim.

c = calcSoundSpeed(22);
L = size(ts.data,1);            % number of samples in pulse
M = size(sLoc,1);               % number of sensor channels

% calculate range and time delay from source to each sensor
range = dist(sLoc',source');    % Euclidean range to source
delay = range./c;               % time delay to source

t0 = ceil(delay * ts.fs);       % delay to source in samples
TMax = max(t0)+L;               % maximum time period for simulation
%ceil(max(delay) + L);

% init frame matrix
sim.fs = ts.fs;
sim.time = (0:TMax-1)/sim.fs;
sim.data = zeros(TMax,M);

% add signal to each channel with appropriate delay
for m=1:M
    sim.data(t0(m):t0(m)+L-1,m) = ts.data;
    
    % delay with group delay FIR filter
    %b = [zeros(1,t0(m)-1) 1];
    %sim.data(:,m) = fftfilt(b,[ts.data; zeros(TMax-L,1)]);
    
    % delay by phase shifting signal
    %sim.data(:,m) =  ts.data * exp(-j*2*pi);
end
