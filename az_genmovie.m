% script to generate movie of single beam over frequency

beamfile = 'starbuck_beams_01.mat';
prefix = 'starbuck';

load(beamfile)

for i=1:numel(beam)
    
    vidname = [prefix '_' sprintf('%.3d',callIdx(i)) '.avi'];
    
    % init AVI file
    %aviobj = avifile(vidname);
    %aviobj.Quality = 100;
    %aviobj.Compression = 'None';
    %aviobj.Fps = 15;
    vidobj = VideoWriter(vidname);
    vidobj.FrameRate = 10;
    
    open(vidobj);
    
    % open new figure
    fh = figure('MenuBar','none','ToolBar','none','NumberTitle','off');
    
    % iterate over frequency
    for j=1:numel(beam{i}.f)
        plotBeamPattern(beam{i},beam{i}.f(j),'surf','fft',fh);
        
%        pause
        drawnow
        
        % capture and add frame
        F = getframe(fh);
        %aviobj = addframe(aviobj, F(j));
        writeVideo(vidobj,F);
        
    end
    
    % close movie file
    close(fh)
    %aviobj = close(aviobj);
    close(vidobj);
    
    
    
end
