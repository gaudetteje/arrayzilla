function plotTimeSeries(fname1,fname2,callmap,varargin)
%  PLOTTIMESERIES  plots time series data from entire recording

% optional inputs
call_index = [];
if nargin > 3
    call_index = varargin{1};
end

% read all data from selected channels
data1 = read_SRZ(fname1,[],[72 1 52]) * 5 * 2^-16;
data2 = read_SRZ(fname2,[],[54 1]) * 5 * 2^-16;


% init indices
idx = [];
t = []; %zeros(1,);
t0 = zeros(1,length(callmap));
t1 = zeros(1,length(callmap));

% iterate over each call to build sample/time indices
for k = 1:length(callmap)
    c = callmap(k);
    
    % append to sample index
    idx = [idx (c.s0:c.s1)];
    
    % append to time index
    t = [t linspace(c.t0, c.t1, (c.s1-c.s0)+1)];
    
    % construct array of markers
    t0(k) = c.t0;
    t1(k) = c.t1;
    
end

data = [data1(idx,:) data2(idx,:)];
clear data1 data2 idx

figure
plot(t, data)
grid on; hold on;
legend('Center','Top Left','Top Right','Bottom Left','Bottom Right')
plot(t0,zeros(size(t0)),'.b')
plot(t1,zeros(size(t1)),'.r')
xlabel('Time (sec)')
ylabel('Amplitude (Volts)')

if ~isempty(call_index)
    plot(callmap(call_index).t0, zeros(size(call_index)), 'bo')
    plot(callmap(call_index).t1, zeros(size(call_index)), 'ro')
end
