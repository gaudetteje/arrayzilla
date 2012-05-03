function beam = az_calcbeam(fd,array,src,varargin)
% AZ_CALCBEAM  calculates and interpolates beams across frequency

fprintf('\n\n*****************************************\n')
fprintf('Computing beam pattern\n\n')

% optional parameters
DEBUG = 0;
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
fRng = (15:10:95)*1e3;     % define frequency range bins
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


%% interpolate spectral density (FFT) data points

% iterate over each frequency range
for i = 1:length(fRng)-1

    % consider only frequency samples within range
    idx = sparse(fd.freq >= fRng(i) & fd.freq < fRng(i+1));

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

%% interpolate Hilbert Spectral Analysis (MCA) data points
if isfield(fd,'if1')
    clear res
    
    % iterate over each frequency range
    for i = 1:length(fRng)-1
        
        % consider only frequency samples within range
        idx = sparse(fd.if1 >= fRng(i) & fd.if1 < fRng(i+1));
        
        % extract magnitudes
        for ch = 1:size(idx,2)
            res(ch) = db(mean(fd.ia1(idx(:,ch))));
        end
        
        % interpolate data on uniform grid
        if exist('TriScatteredInterp','file')
            B = TriScatteredInterp(xPos(:),yPos(:),res(:),INTERPMODE);
            beam.HSA1(:,:,i) = B(beam.X,beam.Y);
        else
            warning('TriScatteredInterp is missing!  Using griddata instead')
            beam.HSA1(:,:,i) = griddata(beam.X, beam.Y, res, xPos, yPos, 'cubic');
        end
        
        if DEBUG
            fh = figure;
            surf(beam.X,beam.Y,beam.HSA1(:,:,i))
            pause, 
            drawnow, close(fh)
        end
    end
    
end

if isfield(fd,'if2')
    clear res
    
    % iterate over each frequency range
    for i = 1:length(fRng)-1
        
        % consider only frequency samples within range
        idx = sparse(fd.if2 >= fRng(i) & fd.if2 < fRng(i+1));
        
        % extract magnitudes
        for ch = 1:size(idx,2)
            res(ch) = db(mean(fd.ia2(idx(:,ch))));
        end
        
        % interpolate data on uniform grid
        if exist('TriScatteredInterp','file')
            B = TriScatteredInterp(xPos(:),yPos(:),res(:),INTERPMODE);
            beam.HSA2(:,:,i) = B(beam.X,beam.Y);
        else
            warning('TriScatteredInterp is missing!  Using griddata instead')
            beam.HSA2(:,:,i) = griddata(beam.X, beam.Y, res, xPos, yPos, 'cubic');
        end
        
        if DEBUG
            fh = figure;
            surf(beam.X,beam.Y,beam.HSA1(:,:,i))
            pause, 
            drawnow, close(fh)
        end
    end
    
end
