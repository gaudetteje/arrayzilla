function plotTimeSeries(fname1,fname2,callmap)
% PLOTTIMESERIES  plots time series data from entire recording
%
% plotTimeSeries(FNAME1,FNAME2,CALLMAP) plots the raw time series data for
% each call contained in CALLMAP


N = numel(callmap);            % total number of events
ch1 = [72 1 110];
ch2 = [49 1];

% init indices
data1 = cell(N,1);
data2 = cell(N,1);
t1 = cell(N,1);
t2 = cell(N,1);

figure
hold on;

% iterate over each call
for k = 1:N
    c = callmap(k);
    
    % sample index
    idx1 = (c.s0(1) : c.s1(1));
    idx2 = (c.s0(2) : c.s1(2));
    
    % read all data from selected channels
    data1{k} = read_SRZ(fname1,idx1,ch1) * (5 * 2^-16);
    data2{k} = read_SRZ(fname2,idx2,ch2) * (5 * 2^-16);
    
    % append to time index
    t1{k} = linspace(c.t0(1), c.t1(1), (c.s1(1)-c.s0(1))+1);
    t2{k} = linspace(c.t0(2), c.t1(2), (c.s1(2)-c.s0(2))+1);
    
    % plot current call
    plot(t1{k}, data1{k})
    plot(t2{k}, data2{k})

end

grid on
xlabel('Time (sec)')
ylabel('Amplitude (Volts)')
