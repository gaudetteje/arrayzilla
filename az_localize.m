function src = az_localize(a, ts, varargin)
% AZ_LOCALIZE  localize sound source for a single pulse and map array
% coordinates to angular map (az/el)
%
% SRC = az_localize(ARRAY,TS) uses the time series data, TS, to localize
%   a point source given the array coordinates in the ARRAY struct
%
% SRC = az_localize(ARRAY,TS,true) plots the array with original and
%   projected angular coordinates about the source origin
%
% SRC = az_localize(ARRAY,TS,false,false) bypasses the localization and
%   hard codes the source position for debugging purposes
%
% Note:  The (0,0,Z) vector is normal to the array plane assuming the point
%   source position is the origin.

fprintf('\n***********************************************\n')
fprintf('Localizing sound sources and calculating angles\n')

% flags
DEBUG = false;      % generates simulated data for debugging
LOCMODE = true;     % attempts localization of sound source using TDOA
VERBOSE = false;    % displays lengthy information about beam resolution/swath

% optional inputs
PLOTFLAG = false;
switch nargin
    case 2
    case 3
        PLOTFLAG = varargin{1};
    case 4
        PLOTFLAG = varargin{1};
        LOCMODE = varargin{2};
    otherwise
        error('Incorrect number of parameters entered')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beamform data using TDOA method


% select subset of microphones for TDOA
ch = [110 83 28 56 58 85  1 105 106 79 33 62 63 100 101 94 10 11 12 13 70 14 43 44 33 34  7 78 107 108 53  1 85 58 83 28];  %74
bd = [  1  1  1  1  2  2  2   1   1  1  2  2  2   1   1  2  2  2  1  1  2  2  2  2  1  1  1  2   2   2  2  1  1  1  2  2];  %1
[idx,bad] = az_chanindex(ch, bd, a);
idx(bad) = [];          % remove missing or known bad channels
micpos = [a.xPos(idx); a.yPos(idx); zeros(1,length(idx))];  %1e-3*randn(1,length(idx))];         % use design coordinates
%load('arrayzilla_coords.mat','micpos')                                 % use sectioned coordinates
%micpos(:,bad) = [];     % remove missing or known bad channels

if DEBUG
%%%%%%%%%%%%%%%%
    % sim parameters
    fs = 236660.3;                   % sampling rate of generated waveform [Hz]
    srcLoc =  [.9 .7 1.2];      % (x,y,z) coordinates

    % create time series of emitted pulse
    f0 = 105e3;             % initial frequency [Hz]
    f1 = 15e3;              % final frequency [Hz]
    T = 0.001;              % pulse length [sec]
    N0 = 10^(-Inf/20);          % add uncorrelated noise

    t = (0:1/fs:T)';
    IF = f0*ones(size(t)) + (f1-f0)/(T).*t;
    IA = raisedcos(length(IF));
    sim.time = t;
    sim.data = real(gen_ifpulse(fs,IF,0,IA));
    sim.fs = fs;

    sim=array_sim(sim, micpos', srcLoc);

    coords = TDOA_frame(sim.data', micpos', sim.fs );

else
    if LOCMODE
        coords = TDOA_frame(ts.data(:,idx)', micpos', ts.fs);
    else
        warning('AZ_LOCALIZE:tdoa','Hard coding position due to inconsistent results');
        coords = [0.93 0.78 1 0];
    end
end

src.xSrc = coords(1);
src.ySrc = coords(2);
src.zSrc = coords(3);
src.residual = coords(4:end);


fprintf('\nEstimated source location is (%g, %g, %g) [m]\n', src.xSrc, src.ySrc, src.zSrc)

% redefine array points around source origin
xPos = a.xPos - src.xSrc;
yPos = a.yPos - src.ySrc;
zPos = ones(size(a.xPos)) * src.zSrc;

% calculate source to mic distance and angle for each channel
[src.az, src.el, src.rng] = cart2sph(xPos, zPos, yPos);
src.az = 90 - src.az*180/pi;        % rotate by -90 and convert to degrees
src.el = src.el*180/pi;             % convert to degrees


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% print calculated array parameters
if VERBOSE
    azRes = diff(src.az); azRes(azRes <= 0) = [];
    elRes = diff(src.el); elRes(elRes <= 0) = [];

    fprintf('\nBeam Coverage:\n')
    fprintf('Min. / Max. Horizontal Angle\n')
    fprintf('   %g / %g degrees\n', [min(src.az) max(src.az)]);
    fprintf('Min. / Max. Vertical Angle\n')
    fprintf('   %g / %g degrees\n\n', [min(src.el) max(src.el)]);

    fprintf('Beam Resolution:\n')
    fprintf('Min. / Max. Horizontal Resolution\n')
    fprintf('   %g / %g degrees\n', [min(azRes) max(azRes)]);
    fprintf('Min. / Max. Vertical Resolution\n')
    fprintf('   %g / %g degrees\n\n', [min(elRes) max(elRes)]);

    fprintf('Minimum / Maximum Euclidean Distance\n')
    fprintf('  %g / %g meters\n', [min(src.rng) max(src.rng)]);
end

% plot spatial representation with source locations (Note: Y and Z are reversed to simplify 3d rotation)
if PLOTFLAG
    %% generate plot of planar array

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

    
    
%     figure
%     plot3(micpos(1,:),micpos(3,:),micpos(2,:),'.')      % sensor locations
%     grid on
%     hold on
%     set(gca,'YDir','reverse')
%     xlabel('x [m]')
%     ylabel('z [m]')
%     zlabel('y [m]')
%     plot3(xSrc,zSrc,ySrc,'r.')           % estimated location
%     title(sprintf('Source @ (%.3f,%.3f,%.3f)',xSrc,ySrc,abs(zSrc)));
%     axis equal
    
    drawnow
end
