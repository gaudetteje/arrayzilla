function genVideoTrack(beam,cdata,avifile,varargin)
% GENVIDEOTRACK  compiles beam pattern images into an AVI file
%
% genVideoTrack(BEAMS,CALLDATA,FILENAME) takes beam data in an Nx1 cell
%    array and reconstructs the beam progression over time at a constant
%    frame rate.
%
% NOTE:  MATLAB has a known bug in Windows Vista (and 7) with the OpenGL
% renderer.  We can use zbuffer or run software graphics processing.  The
% latter will be slower, but zbuffer doesn't support transparency.


close all

if ispc
    opengl('software');
end

% default parameters
D = 5;                          % slow audio playback by factor of D
tBuf = 0.1;                     % add short buffer to beginning
fRate = 30;                     % specify frame rate
f = 60e3;                       % frequency bin to plot

% override frequency parameter
if nargin > 3
    f = varargin{1};
end

% compute timescale and frame numbers
t = [cdata(:).t0];
t = tBuf + t -  min(t);         % call times (original rate)
T = tBuf + (max([cdata.t1]) - min([cdata.t0]));        % total duration

t2 = t * D;                     % relative call times (slowed by D with padding)
T2 = T * D;                     % call times (playback rate)
c = floor(t2*fRate);            % call times (frame number)
N = ceil(T2*fRate);             % total number of frames



% initialize AVI video
vidobj = VideoWriter(avifile);
vidobj.FrameRate = fRate;
open(vidobj);

% open new figure
hFig = figure('MenuBar','none','ToolBar','none','NumberTitle','off','NextPlot','replace');
set(gcf, 'Renderer', 'zbuffer')

% plot transparent placeholder
cNum = 1;
cNext = c(cNum);
c(end+1) = N+1;                 % append after last call to avoid error
hAxis = plotBeamPattern(beam{cNum},f,'surf','fft',hFig);
hSurf = findobj(hAxis,'type','surface');
set(hSurf,'edgealpha',0);
set(hSurf,'facealpha',0);
hCont = findobj(hAxis,'type','patch');
set(hCont,'edgealpha',0);

% iterate over each video frame
for n=1:N
    
    % determine if we are plotting a new beam or updating a previous beam
    if (n >= cNext)
        
        % plot new beam on existing figure
        if ~isempty(beam{cNum})     % skip beam if empty
            hAxis = plotBeamPattern(beam{cNum},f,'surf','fft',hFig);
        %else
            % just plot previous beam for now - NEED TO FIX FAILED BEAM
            % PROCESSING!
            %hAxis = plotBeamPattern(beam{cNum-1},f,'surf','fft',hFig);
            %alpha('clear') % remove previous beam
        end
        
        % get surface handle from axis
        hSurf = findobj(hAxis,'type','surface');
        
        % update call number
        cNum = cNum + 1;
        cNext = c(cNum);
        
    else

        % reduce transparency by percentage to give exponential fading
        a = get(hSurf,'FaceAlpha');
        set(hSurf,'FaceAlpha',a*0.95);
        
    end
    
    title(sprintf('%.2f sec., Call %d, %g kHz', n/(D*fRate), cNum-1, f*1e-3))
    
    % capture and add frame
    refreshdata
    drawnow
    
    F = getframe(hFig);
    writeVideo(vidobj,F);
    
end

% close movie file
close(hFig)
close(vidobj);
