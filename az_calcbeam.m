function beam = az_calcbeam(fd,array,src,varargin)
% AZ_CALCBEAM  calculates and interpolates beams across frequency

fprintf('\n\n*****************************************\n')
fprintf('Computing beam pattern\n\n')

% optional parameters
AXISMODE = 'ang';  %'pos'
INTERPMODE = 'natural'; %'nearest','linear';
if nargin > 3
    AXISMODE = varargin{1};
    if nargin > 4
        INTERPMODE = varargin{2};
    end
end


%% hard coded parameters

% frequency & grid spacing
fRng = (35:10:95)*1e3;     % define frequency range bins
%fRng = (12.5:5:102.5)*1e3;    % define frequency range bins
xPts = 19*3;            % numebr of x grid points
yPts = 12*3;            % number of y grid points


%% assign (x,y) data grid indices
switch AXISMODE(1:3)
    case 'ang'
        xPos = src.az;
        yPos = src.el;
        xStr = 'Azimuth (deg)';
        yStr = 'Elevation (deg)';
    case 'pos'
        xPos = array.xPos;
        yPos = array.yPos;
        xStr = 'Array Location (m)';
        yStr = 'Array Location (m)';
    otherwise
        error('Unknown  mode')
end

% define tightly sampled uniform grid
beam.f = fRng(1:end-1) + diff(fRng)/2;
beam.x = linspace(min(xPos), max(xPos), xPts);
beam.y = linspace(min(yPos), max(yPos), yPts);
[beam.X, beam.Y] = meshgrid(beam.x,beam.y);


%% iterate over each frequency range
for i = 1:length(fRng)-1
    
    % consider only frequency samples within range
    idx = fd.freq >= fRng(i) & fd.freq < fRng(i+1);
    
    % average the data points
    res = mean(fd.magdb(idx,:),1);
    
    % interpolate data on uniform grid
    if exist('TriScatteredInterp','file')
        B = TriScatteredInterp(xPos(:),yPos(:),res(:),INTERPMODE);
        beam.Z(:,:,i) = B(beam.X,beam.Y);
    else
        warning('TriScatteredInterp is missing!  Using griddata instead')
        beam.Z(:,:,i) = griddata(beam.X, beam.Y, res, xPos, yPos, 'cubic');
    end
    
end
