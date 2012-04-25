function a = az_positions(varargin)
% AZ_POSITIONS  defines a uniform grid of array sensor positions and dimensions
%
% ARRAY = AZ_POSITIONS([Nx Ny],[dx dy]) returns a struct containing the array
%     sensor positions in a Nx by Ny matrix with pitch (dx,dy) in meters
% ARRAY = AZ_POSITIONS([Nx Ny],[dx dy],[x0 y0]) recenters the array at
%     the position index (x0,y0)
% ARRAY = AZ_POSITIONS  assumes the default set of parameters below
%
% Inputs:
%     [Nx Ny] = Number of sensor [int]
%     [dx dy] = Pitch of sensor spacing in meters [double]
%     [x0 y0] = Origin index of first sensor [int]
%
% Defaults:
%     Nx = 3;
%     Ny = 3;
%     dx = 0.01;
%     dy = 0.01;
%     x0 = (Nx+1)/2;
%     y0 = (Ny+1)/2;
%     
% Example:
%    array = az_positions([19 12],[4 5]*.0254,[.5 .5])
%

% optional parameters
PLOTFLAG = false;

% get optional input parameters
switch nargin
    case 0
        Nx = 3;
        Ny = 3;
        dx = 0.01;
        dy = 0.01;
        x0 = 2;
        y0 = 2;
    case 1
        Nx = varargin{1}(1);
        Ny = varargin{1}(2);
        dx = 0.01;
        dy = 0.01;
        x0 = (Nx+1)/2;
        y0 = (Ny+1)/2;
    case 2
        Nx = varargin{1}(1);
        Ny = varargin{1}(2);
        dx = varargin{2}(1);
        dy = varargin{2}(2);
        x0 = (Nx+1)/2;
        y0 = (Ny+1)/2;
    case 3
        Nx = varargin{1}(1);
        Ny = varargin{1}(2);
        dx = varargin{2}(1);
        dy = varargin{2}(2);
        x0 = varargin{3}(1);
        y0 = varargin{3}(2);
end


%% Calculate uniform mic (X,Y) positions in metric units [m]
a.xPos = repmat(((1:Nx)-x0)*dx,1,Ny);
a.yPos = sort(repmat(((1:Ny)-y0)*dy,1,Nx));

% plot results
if PLOTFLAG
    figure
    plot(a.xPos,a.yPos,'r.');
    title('Microphone positions for Arrayzilla')
    xlabel('X position (m)')
    ylabel('Y position (m)')
    grid on
    hold on
    
    plot(0,0,'k+')            % origin of absolute array
    axis equal
    axis([min(a.xPos - dx/2) max(a.xPos + dx/2) ...
          min(a.yPos - dy/2) max(a.yPos + dy/2)])
    
    drawnow
end
