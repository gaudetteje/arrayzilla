function plotTriggerEvents(hdr,events)
% PLOTTRIGGEREVENTS  plots auxiliary channel data vs. time
%
% plotTriggerEvents(HDR,EVENTS) takes all events contained in EVENTS and
% plots the raw auxiliary "sync" channel data saved to the HDR struct

N = numel(events);            % total number of events

% init indices
data1 = cell(N,1);
data2 = cell(N,1);
t1 = cell(N,1);
t2 = cell(N,1);

figure
hold on;

% iterate over each call
for k = 1:N
    c = events(k);
    
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
