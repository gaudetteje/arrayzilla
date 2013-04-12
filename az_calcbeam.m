function beam = az_calcbeam(fd,src,array,varargin)
% AZ_CALCBEAM  calculates and interpolates beams across frequency

fprintf('\n***********************************************\n')
fprintf('Interpolating and smoothing beam pattern\n')

% optional parameters
DEBUG = 0;
INTERPMODE = 'natural'; %'nearest','linear';
if nargin > 3
    INTERPMODE = varargin{1};
end



%% hard coded parameters

% frequency & grid spacing
fRng = (9.5:1:110.5)*1e3;      % define frequency range bins
azPts = 19*5;                   % number of x grid points
elPts = 12*5;                   % number of y grid points


%% assign (az,el) data grid indices
az = src.az;
el = src.el;

% define tightly sampled uniform grid
beam.f = fRng(1:end-1) + diff(fRng)/2;
beam.az = linspace(min(az), max(az), azPts);
beam.el = linspace(min(el), max(el), elPts);
[beam.AZ, beam.EL] = meshgrid(beam.az,beam.el);


%% interpolate spectral density (FFT) data points

% iterate over each frequency range
for i = 1:length(fRng)-1

    % consider only frequency samples within range
    idx = sparse(fd.freq >= fRng(i) & fd.freq < fRng(i+1));

    % average the data points
    res = mean(fd.magdb(idx,:),1);

    % interpolate data on uniform grid
    if exist('TriScatteredInterp','file')
        B = TriScatteredInterp(az(:),el(:),res(:),INTERPMODE);
        beam.FFT(:,:,i) = B(beam.AZ,beam.EL);
    else
        warning('TriScatteredInterp is missing!  Using griddata instead')
        beam.FFT(:,:,i) = griddata(beam.AZ, beam.EL, res, az, el, 'cubic');
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
            B = TriScatteredInterp(az(:),el(:),res(:),INTERPMODE);
            beam.HSA1(:,:,i) = B(beam.AZ,beam.EL);
        else
            warning('TriScatteredInterp is missing!  Using griddata instead')
            beam.HSA1(:,:,i) = griddata(beam.AZ, beam.EL, res, az, el, 'cubic');
        end
        
        if DEBUG
            fh = figure;
            surf(beam.AZ,beam.EL,beam.HSA1(:,:,i))
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
            B = TriScatteredInterp(az(:),el(:),res(:),INTERPMODE);
            beam.HSA2(:,:,i) = B(beam.AZ,beam.EL);
        else
            warning('TriScatteredInterp is missing!  Using griddata instead')
            beam.HSA2(:,:,i) = griddata(beam.AZ, beam.EL, res, az, el, 'cubic');
        end
        
        if DEBUG
            fh = figure;
            surf(beam.AZ,beam.EL,beam.HSA1(:,:,i))
            pause, 
            drawnow, close(fh)
        end
    end
end
