function plotTriggerEvents(hdr,callmap)
% PLOTTRIGGEREVENTS  plots trigger data from entire recording
%
% plotTriggerEvents(HDR,CALLMAP) takes all calls contained in CALLMAP and
% plots the raw auxiliary "sync" channel data

N = numel(callmap);            % total number of events

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
    
    % plot sync channel
    plot(hdr(1).time(idx1), hdr(1).aux(idx1))
    plot(hdr(2).time(idx2), hdr(2).aux(idx2),'--r')
end

grid on
xlabel('Time (sec)')
ylabel('Amplitude (Volts)')
