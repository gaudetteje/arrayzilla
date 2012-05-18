function plotArrayPoints(a,src)
% PLOTARRAYPOINTS  plots a 3D diagram showing the planar array sampling
%   points and the angular coordinates they map to.
%
%

close all

%% generate plot of planar array

% redefine array points around source origin
xPos = a.xPos - src.xSrc;
yPos = a.yPos - src.ySrc;
zPos = src.zSrc * ones(size(xPos));

% swap y and z for plot rotation
figure
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

%% rotate 
N = 200;
az0 = linspace(-90,-270,N);
el0 = [linspace(90,0,N/2) zeros(1,N/2)];
for n = 1:N
    view(az0(n),el0(n))
    %drawnow
    pause(.001)
end
