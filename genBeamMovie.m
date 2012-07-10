function genBeamMovie(beam,avifile,varargin)
% GENBEAMMOVIE  creates an AVI video showing a single beam over frequency

% force beam back into cell array
if ~iscell(beam)
    beam = {beam};
end


for i=1:numel(beam)
    
    % init AVI file
    vidobj = VideoWriter(avifile);      % need to include call number here, otherwise it will be overwritten
    vidobj.FrameRate = 10;
    
    open(vidobj);
    
    % open new figure
    fh = figure('MenuBar','none','ToolBar','none','NumberTitle','off');
    
    % iterate over frequency
    for j=1:numel(beam{i}.f)
        
        plotBeamPattern(beam{i},beam{i}.f(j),'surf','fft',fh);
        drawnow
        
        % capture and add frame
        F = getframe(fh);
        writeVideo(vidobj,F);
        
    end
    
    % close movie file
    close(fh)
    close(vidobj);
    
end
