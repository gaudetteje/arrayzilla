% script to process data separated by trials

prefix = 'starbuck';
callfile = sprintf('%s_callmap.mat',prefix);
fname1 = sprintf('%s_side1.srz',prefix);
fname2 = sprintf('%s_side2.srz',prefix);

DEBUG = false;

load(callfile)

% extract starting times of all events
t = [callmap(:).t0];

% use difference between events as trial indicator
t = find(diff(t(1:2:end)) > 0.6) + 2;

% remove trials with less than 50 pulses
t = t(diff(t) > 50);

% add first and last events
t = [1 t numel(callmap)+1];


% iterate over each trial
tic
for n = 1:numel(t)-1
    try
        if DEBUG
            plotTimeSeries('starbuck_side1.srz','starbuck_side2.srz',callmap(t(n):t(n+1)-1))
        end

        % process each trial
        beamname = sprintf('%s_beams_%d-%d',prefix,t(n),t(n+1)-1);
        [b, ref] = az_process_beams(fname1,fname2,t(n):t(n+1)-1,beamname);
        
        wavname = sprintf('%s_call_%d-%d.wav',prefix,t(n),t(n+1)-1);
        genAudioTrack(ref, wavname);
        
    catch
        warning('Crap!')
        continue
    end
end
toc

