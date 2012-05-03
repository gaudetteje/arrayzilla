function src = az_localize(ts, a, varargin)
% AZ_LOCALIZE  localize sound sources for each pulse and map array
% coordinates to angular map (az/el)
%
% Note:  (0,0) vector is normal to array plane

fprintf('\n\n***********************************************\n')
fprintf('Localizing sound sources and calculating angles\n')

DEBUG = false;

% optional inputs
PLOTFLAG = false;
if nargin > 2
    PLOTFLAG = varargin{1};
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    coords = TDOA_frame(ts.data(:,idx)', micpos', ts.fs);
end

xSrc = 0.93; %coords(1);
ySrc = 0.78; %coords(2);
zSrc = 1; %coords(3);
warning('Hard coding position due to inconsistent results') %warning('FORCING Z=1m')

fprintf('\nEstimated source location is (%g, %g, %g) [m]\n',xSrc,ySrc,zSrc)

% calculate source to mic angle
az = atan2(a.xPos-xSrc,zSrc)*180/pi;
el = atan2(a.yPos-ySrc,zSrc)*180/pi;

% calculate source to mic distance for each channel
rng = dist([xSrc; ySrc; zSrc], [a.xPos; a.yPos; zeros(size(a.ch))]);


% return struct of results
src.xSrc = xSrc;
src.ySrc = ySrc;
src.zSrc = zSrc;
src.residual = coords(4:end);

src.az = az;
src.el = el;
src.rng = rng;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% print calculated array parameters
azRes = diff(az); azRes(azRes <= 0) = [];
elRes = diff(el); elRes(elRes <= 0) = [];

fprintf('\nBeam Coverage:\n')
fprintf('Min. / Max. Horizontal Angle\n')
fprintf('   %g / %g degrees\n', [min(az) max(az)]);
fprintf('Min. / Max. Vertical Angle\n')
fprintf('   %g / %g degrees\n\n', [min(el) max(el)]);

fprintf('Beam Resolution:\n')
fprintf('Min. / Max. Horizontal Resolution\n')
fprintf('   %g / %g degrees\n', [min(azRes) max(azRes)]);
fprintf('Min. / Max. Vertical Resolution\n')
fprintf('   %g / %g degrees\n\n', [min(elRes) max(elRes)]);

fprintf('Minimum / Maximum Euclidean Distance\n')
fprintf('  %g / %g meters\n', [min(rng) max(rng)]);


% plot spatial representation with source locations (Note: Y and Z are reversed to simplify 3d rotation)
if PLOTFLAG
    figure
    plot3(micpos(1,:),micpos(3,:),micpos(2,:),'.')      % sensor locations
    grid on
    hold on
    set(gca,'YDir','reverse')
    xlabel('x [m]')
    ylabel('z [m]')
    zlabel('y [m]')
    plot3(xSrc,zSrc,ySrc,'r.')           % estimated location
    title(sprintf('Source @ (%.3f,%.3f,%.3f)',xSrc,ySrc,abs(zSrc)));
    axis equal
    
    drawnow
end
