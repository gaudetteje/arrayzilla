function plotAuxChan(varargin)
% PLOTAUXCHAN  plots auxiliary channel data vs. time
%
% plotAuxChan(HDR,EVENTS) takes all events contained in EVENTS and
%    plots the raw auxiliary "sync" channel data saved to the HDR struct
% plotAuxChan(FNAME1,FNAME2,EVENTS) reads header info from each event in
%    the files specified
% 

switch nargin
    case 3
        fname1 = varargin{1};
        fname2 = varargin{2};
        events = varargin{3};
    case 2
        hdr = varargin{1};
        events = varargin{2};
    otherwise
        error('Incorrect number of parameters entered')
end

N = numel(events);              % total number of events
M = max(events(end).s1);        % upper bound on number of events

% read data from file if necessary
if ~exist('hdr','var')
    % prealloc struct
    [hdr(1:2).aux] = deal(zeros(1,M));
    [hdr(1:2).time] = deal(zeros(1,M));
    
    % iterate over each event
    for n = 1:N
        idx1 = (events(n).s0(1) : events(n).s1(1));
        hdr(1).aux(idx1) = read_SRZ_header_field(fname1,idx1,9);
        hdr(1).time(idx1) = linspace(events(n).t0(1), events(n).t1(1), numel(idx1));
        
        idx2 = (events(n).s0(2) : events(n).s1(2));
        hdr(2).aux(idx2) = read_SRZ_header_field(fname2,idx2,9);
        hdr(2).time(idx2) = linspace(events(n).t0(2), events(n).t1(2), numel(idx2));
    end
end

% plot results
figure
hold on;

% iterate over each event
for k = 1:N
    c = events(k);
    
    % sample index
    idx1 = (c.s0(1) : c.s1(1));
    idx2 = (c.s0(2) : c.s1(2));
    
    % plot aux channel
    plot(hdr(1).time(idx1), hdr(1).aux(idx1))
    plot(hdr(2).time(idx2), hdr(2).aux(idx2),'--r')
end

grid on
xlabel('Time (sec)')
ylabel('Amplitude (Volts)')
