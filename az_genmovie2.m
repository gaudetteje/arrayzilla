clear 
close all
clc

% script to generate movie of multiple beams over time

beamfile = 'starbuck_beams_1-325.mat';%325-450.mat';
callfile = 'starbuck_callmap.mat';
prefix = 'starbuck';

load(beamfile);
load(callfile);


%% parameters
f = 68e3;                               % select single frequency to plot
T0 = callmap(callIdx(1)).t0(1);         % get absolute time of first call



% init AVI file
vidname = sprintf('%s_%.1dkHz_calls_%d-%d.avi', prefix, f*1e-3, callIdx(1), callIdx(end));
vidobj = VideoWriter(vidname);
vidobj.FrameRate = 10;
open(vidobj);

% open new figure
fh = figure('MenuBar','none','ToolBar','none','NumberTitle','off');

% iterate over beam number
for i=1:numel(beam)

    % get absolute time from first call
    
    
    if ~isempty(beam{i})
        plotBeamPattern(beam{i},f,'surf','fft',fh);
        title(sprintf('Call %d - %.2fms',callIdx(i), (callmap(callIdx(i)).t0(1) - T0) * 1e3 ))

        %        pause
        drawnow

        % capture and add frame
        F = getframe(fh);
        %aviobj = addframe(aviobj, F(j));
        writeVideo(vidobj,F);
    end
    
end

% close movie file
close(fh)
%aviobj = close(aviobj);
close(vidobj);