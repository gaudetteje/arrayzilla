function plotArrayPoints(a,varargin)
% PLOTARRAYPOINTS  plots a 3D diagram showing the planar array sampling
%   points and the angular coordinates they map to.
%
% plotArrayPoints(ARRAY,SOURCE) plots array geometry and mapped angular
%   coordinates relative to the sound source
%
% plotArrayPoints(ARRAY) can be used when no source information is
%   available.  In this case, the source is assumed to be at 1 meter normal
%   to the center of the array
%
% plotArrayPoints(..,true) will generate a video circumventing the array
%
% plotArrayPoints(..,true,VIDNAME) saves the video to the specified name

close all
GENMOVIE = false;
avifile = 'rotate_array_points.avi';

switch(nargin)
    case 1
    case 2
        if islogical(varargin{1})
            GENMOVIE = varargin{1};
        else
            src = varargin{1};
        end
    case 3
        if islogical(varargin{1})
            GENMOVIE = varargin{1};
            avifile = varargin{2};
        else
            src = varargin{1};
            GENMOVIE = varargin{2};
        end
    case 4
        src = varargin{1};
        GENMOVIE = varargin{2};
        avifile = varargin{3};
    otherwise
        error('Incorrect number of parameters entered')
end

% if no source specified, force origin to 1 meter normal to center of array
if ~exist('src','var')
    src.xSrc = (max(a.xPos)-min(a.xPos))/2;
    src.ySrc = (max(a.yPos)-min(a.yPos))/2;
    src.zSrc = 1;
end

% assume the following spacing and dimensions if not specified
if ~isfield(a,'dx')
    a.Nx = 19;
    a.Ny = 12;
    a.dx = 0.1016;
    a.dy = 0.1270;
    a.x0 = 0.5000;
    a.y0 = 0.5000;
end

%% generate plot of planar array

% redefine array points around source origin
xPos = a.xPos - src.xSrc;
yPos = a.yPos - src.ySrc;
zPos = src.zSrc * ones(size(xPos));

% swap y and z for plot rotation
fh = figure('color','white');
plot3(xPos,zPos,yPos,'r.');
title('Microphone positions for Arrayzilla')
xlabel('X position (m)')
ylabel('Z position (m)')
zlabel('Y position (m)')
grid on
hold on

% set axes
axis equal
view(-134,46)

% define linear array boundaries
bL = min(xPos - .0127);
bR = max(xPos + .0127);
bD = min(yPos - .0127);
bU = max(yPos + .0127);

% plot array boundary
plot3([bL bR bR bL bL], src.zSrc*ones(1,5), [bD bD bU bU bD], ...
    'k', 'linewidth', 2)
axis([bL-a.dx bR+a.dx 0 src.zSrc bD-a.dy bU+a.dy])

plot3(0,src.zSrc,0,'k.')            % mark the array point normal to origin

% set tick marks
dt = 0.2;
set(gca,'XTick',(dt*floor(bL/dt) : dt : dt*ceil(bR/dt)))
set(gca,'ZTick',(dt*floor(bD/dt) : dt : dt*ceil(bU/dt)))
set(gca,'YTick',(0 : dt : dt*ceil(src.zSrc/dt)))


%% generate spherical coordinates
[az,el,rng] = cart2sph(xPos, zPos, yPos);
rho = 0.4;
[x,y,z] = sph2cart(az,el,rho*ones(size(rng)));

plot3(x,y,z,'b.')

% draw origin
plot3(0,0,0,'k+')

% draw spherical grid around origin
[X,Y,Z] = sphere(360/10);
idx1=1:37;              % elevation
idx2=19:37;             % azimuth
ah = surf(rho*X(idx1,idx2),rho*Y(idx1,idx2),rho*Z(idx1,idx2));
set(ah,'FaceColor','none')
set(ah,'EdgeColor',[.75 .75 .75])

drawnow

%% rotate view and create movie
if GENMOVIE
    % init AVI file
    vidobj = VideoWriter(avifile);      % need to include call number here, otherwise it will be overwritten
    vidobj.FrameRate = 30;
    vidobj.Quality = 25;
    
    open(vidobj);

    % iterate over each view
    N = 200;
    az0 = linspace(-90,-270,N);
    el0 = [linspace(90,0,N/2) zeros(1,N/2)];
    for n = 1:N
        view(az0(n),el0(n))
        pause(.001)
        if GENMOVIE
            % capture and add frame
            drawnow
            F = getframe(fh);
            writeVideo(vidobj,F);
        end
    end

    % clean up
    close(fh)
    close(vidobj)
end