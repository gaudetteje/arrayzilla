function h = plotBeamPattern(B,varargin)
% PLOTBEAMPATTERN  plots an individual measured 3D beam
%
% plotBeamPattern(BEAM) plots the beam pattern across all frequencies
% 
% plotBeamPattern(BEAM,FREQ) plots the beam patten at the specified
%   frequency(s)
%
% plotBeamPattern(BEAM,FREQ,SRC) plots sensor locations over beam data to
%   show where beam data is being interpolated
% 
% plotBeamPattern(...,PLOTMODE) optionally plots beams using a
%   different method (see actual code for available modes)
%
% plotBeamPattern(...,PLOTMODE,DATAMODE) optionally plots beams using a
%   different data set ('fft', 'hsa1', or 'hsa2')
%
% plotBeamPattern(...,PLOTMODE,DATAMODE,FH) forces use of FH as the figure
%   handle
%
% FH = plotBeamPattern(...) also returns the figure handle, FH
%
% BEAM is a struct containing the 3-dimensional array of complex beam data,
% and indices for each dimension.
%
% beam.FFT - i x j x k complex matrix azimuth, elevation, and frequency
% beam.HSA1 - i x j x k complex matrix azimuth, elevation, and frequency
% beam.HSA2 - i x j x k complex matrix azimuth, elevation, and frequency
% beam.el - i-element vector of elevation indices
% beam.az - j-element vector of azimuth indices
% beam.f - k-element vector of frequency indices
% beam.EL - i x j matrix of elevation for each element in Z
% beam.AZ - i x j matrix of azimuth for each element in Z
% 
% To generate X and Y use:
% [beam.AZ, beam.EL] = meshgrid(beam.az, beam.el);


%% set default parameters

% data "massaging" parameters
smMETH = 'box';     % kernel type for smooth3.m
smSIZE = 1;         % kernel size for smooth3.m (set to 1 for disable)

% volume plotting options
%TBD

% surface plotting options
%cMap = 'hot';     %jet;    % colormap
dBrange = 35;           % colorscale depth
zLim = [-100 -20];        % visual depth
dBnorm = true;          % normalize to peak?
azView = 17;            % azimuth angle
elView = 54;            % elevation angle

% contour plotting options
contourLev = -3;        % dB contour level for each frequency line [dB]
colors = {'k','b','g','m','c','r'};


%% optional parameters
FREQ = [];
SRC = [];
PLOTMODE = 'surf'; %'cont'
DATAMODE = 'fft'; %'hsa1'; %
FH = nan;
switch (nargin)
    case 1
    case 2
        FREQ = varargin{1};
    case 3
        FREQ = varargin{1};
        if ischar(varargin{2})
            PLOTMODE = varargin{2};
        else
            SRC = varargin{2};
        end
    case 4
        FREQ = varargin{1};
        if ischar(varargin{2})
            PLOTMODE = varargin{2};
            DATAMODE = varargin{3};
        else
            SRC = varargin{2};
            PLOTMODE = varargin{3};
        end
    case 5
        FREQ = varargin{1};
        SRC = varargin{2};
        PLOTMODE = varargin{3};
        DATAMODE = varargin{4};
    case 6
        FREQ = varargin{1};
        SRC = varargin{2};
        PLOTMODE = varargin{3};
        DATAMODE = varargin{4};
        FH = varargin{5};
    otherwise
        error('Incorrect number of input arguments')
end

% repeat function call over multiple beams, if cell array of beam structs entered
if iscell(B)
    for bNum = 1:numel(B)
        plotBeamPattern(B{bNum},FREQ,SRC(bNum),PLOTMODE,DATAMODE,FH);
    end
    return
end    

fprintf('\n*****************************************\n')
fprintf('Plotting beam pattern\n')

% use specified data set
if strcmp(DATAMODE,'fft')
    B.Z = B.FFT;
elseif strcmp(DATAMODE,'hsa1')
    B.Z = B.HSA1;
elseif strcmp(DATAMODE,'hsa2')
    B.Z = B.HSA2;
else
    error('Unknown DATAMODE parameter.  Should be one of ''fft'', ''hsa1'', or ''hsa2''.')
end

% get indices of desired frequency points
if isempty(FREQ)
    fIdx = 1:numel(B.f);
else
    fIdx = zeros(1,numel(FREQ));
    for i = 1:numel(FREQ)
        res = find(B.f >= FREQ(i), 1, 'first');
        if ~isempty(res)
            fIdx(i) = res;
        end
    end
    fIdx = unique(fIdx);    % remove dups
    fIdx(fIdx == 0) = [];   % remove zeros
end

% set colormap to change with frequency
cBase = gray(64);
colors = jet(numel(B.f));
for i = 1:numel(fIdx)
    cMap{fIdx(i)} = cBase * diag(colors(fIdx(i),:));
end

% figure handle
if ishandle(FH)
    fh = FH;            % get users figure handle
    %if numel(fh) ~= numel(fIdx)
        %fh(end+1) = repmat(fh(end),1,numel(fIdx)-numel(fh));
    %end
else
    fh = nan;%(1,numel(fIdx));    % create array of figure handles to use
end


%% smooth out data before plotting
B.Z = smoothn(B.Z,1e-4);                % interpolate NaN values (especially at missing corners)
B.Z = smooth3(B.Z,smMETH,smSIZE);       % smooth data in all 3 dimensions


%% plot resulting beam patterns
switch PLOTMODE
    
    % plot volumetric plot
    case 'vol'
        
        % convert to volumetric data
        %%% use TriScatteredInterp with 3D uniform meshgrid
%         N = 10;    % number of intensity levels
%         
%         x = repmat(B.AZ,[1 1 N]);
%         y = repmat(B.EL,[1 1 N]);
%         z = zeros(numel(B.el), numel(B.az), N);
%         for i = 1:numel(B.f)
%             z(:,:,i) = B.f(i);
%         end
        
        
        % construct 3D patches for each surface
        for i = fIdx
            if ishandle(fh)
                figure(fh);
            else
                fh = figure;
            end
            fvc = surf2patch(B.AZ,B.EL,B.FFT(:,:,i),B.FFT(:,:,i));
            p{i} = patch(fvc);
            shading faceted;
            set(p{i},'FaceColor','r')
            
            view(3)
            drawnow
            
        end
        
        % see http://www.mathworks.com/matlabcentral/newsreader/view_thread/169205
        
    case 'sph'
        
        % iterate over each frequency beam
        for i = fIdx
            
            % need to convert radius to positive definite units!
            
            % convert (az,el,rho) data points to cartesian coordinates
            [x,y,z] = sph2cart(B.AZ*pi/180, B.EL*pi/180, ones(size(B.EL)));
            %[x,y,z] = sph2cart(B.AZ*pi/180, B.EL*pi/180, B.Z(:,:,i));
%            z = -B.Z(:,:,i) .* cos(B.AZ*pi/180) .* cos(B.EL*pi/180);
%            x = -B.Z(:,:,i) .* sin(B.AZ*pi/180) .* cos(B.EL*pi/180);
%            y = -B.Z(:,:,i) .* sin(B.EL*pi/180);
            
            % find peak value
            dBpeak = max(max(B.Z(:,:,i)));
            fprintf('Peak @ %g kHz = %2.1f dB\n', B.f(i)*1e-3, dBpeak)
            
            % plot interpolated surface
            if ishandle(fh)
                figure(fh);
            else
                fh = figure;
            end
            surf(x,y,z+dBpeak);
%             view(azView,elView)
            shading interp
            lh = light;
            lighting phong
            %lightangle(45,45)
            %axis equal
            %set(gca,'ZLim',[dBpeak-55 dBpeak]);
            
            colormap(cMap);
            if dBnorm
                cRange = [-dBrange 0];
            else
                cRange = [dBpeak-dBrange dBpeak];
            end
            set(gca,'CLim',cRange);
            colorbar
            
            title(sprintf('%g kHz', B.f(i)*1e-3),'fontsize',16)
            xlabel('azim. (\circ)','fontsize',16)
            ylabel('elev. (\circ)','fontsize',16)
            zlabel('magnitude (dB)','fontsize',16)
            
            set(gca,'fontsize',16);
            set(gcf,'color','w');
            
        end
        
        
    % plot 3D surface plots for each frequency bin
    case 'surf'
        
        % iterate over each frequency beam
        for i = fIdx
            
            % find peak value
            dBpeak = max(max(B.Z(:,:,i)));
            fprintf('Peak @ %g kHz = %2.1f dB\n', B.f(i)*1e-3, dBpeak)
            
            % plot interpolated surface
            if ishandle(fh)
                figure(fh);
            else
                fh = figure;
            end
            if dBnorm
                surfc(B.AZ, B.EL, B.Z(:,:,i)-dBpeak);
            else
                surfc(B.AZ, B.EL, B.Z(:,:,i));
            end
            
            % configure view and lighting
            view(azView,elView)
            shading interp
            lh = light;
            lighting phong
            %lightangle(45,45)
            
            % apply transparent mesh grid over surface
            hSurf = findobj('type','surface');
            set(hSurf,'EdgeColor',[.5 .5 .5]);
            set(hSurf,'EdgeAlpha',0.2);
            
            % set colormap
            colormap(cMap{i});
            if dBnorm
                cRange = [-dBrange 0];
            else
                cRange = [dBpeak-dBrange dBpeak];
            end
            set(gca,'CLim',cRange);
            colorbar
            
            % set axes and labels
            if dBnorm
                set(gca,'Zlim',zLim);%+dBpeak);%[dBpeak-50 dBpeak]);
            else
                set(gca,'Zlim',zLim)
            end
            set(gca,'Xlim',[-45 45])
            set(gca,'Ylim',[-45 45])
            
            % adjust contour height
            hCont = findobj('type','patch');
            zd = get(hCont,'ZData');
            for j = 1:numel(hCont)
                set(hCont(j),'ZData',zLim(1)*ones(length(zd{j}),1))
                set(hCont(j),'FaceColor','flat')
            end
            
            
            set(gca,'fontsize',16);
            set(gcf,'color','w');
            
            xlabel('azim. (\circ)','fontsize',16)
            ylabel('elev. (\circ)','fontsize',16)
            zlabel('magnitude (dB)','fontsize',16)
            title(sprintf('%g kHz', B.f(i)*1e-3),'fontsize',16)
            
            
            %pause(0.25)
        end
        
    % plot contours over all frequencies
    case 'cont'
        if ishandle(fh)
            figure(fh);
        else
            fh = figure;
        end
        hold on;
        dBlevels = contourLev * [1 1];
        for i = fIdx
            dBpeak = max(max(B.Z(:,:,i)));
            fprintf('Peak @ %g kHz = %2.1f dB\n', B.f(i)*1e-3, dBpeak)
            contour(B.AZ, B.EL, B.Z(:,:,i),dBpeak+dBlevels,colors{mod(i-1,length(colors))+1},'linewidth',2);
        end
        grid on;
        legend(num2str(1e-3*B.f(fIdx)'))
        
    % plot horizontal slice
    case 'horz'
        k = 28; % select an elevation
        if ishandle(fh)
            figure(fh);
        else
            fh = figure;
        end
        hold on;
        for i = fIdx
            % find peak value
            dBpeak = max(max(B.Z(:,:,i)));
            fprintf('Peak @ %g kHz = %2.1f dB\n', B.f(i)*1e-3, dBpeak)
            
            if dBnorm, zMod = dBpeak; else, zMod = 0; end
            
            plot(B.az, B.Z(k,:,i) - zMod, ...
                colors{mod(i-1,length(colors))+1}, ...
                'linewidth',2)
        end
        grid on;
        
        legend(num2str(1e-3*B.f(fIdx)'))
        xlabel('Azimuth (deg)')
        ylabel('Magnitude (dB)')
        title(sprintf('Horizontal Beam Pattern @ %.1f degrees',B.el(k)))
        
    % plot vertical slice
    case 'vert'
        k = 16; % select an azimuth
        if ishandle(fh)
            figure(fh);
        else
            fh = figure;
        end
        hold on;
        for i = fIdx
            % find peak value
            dBpeak = max(max(B.Z(:,:,i)));
            fprintf('Peak @ %g kHz = %2.1f dB\n', B.f(i)*1e-3, dBpeak)
            
            if dBnorm, zMod = dBpeak; else, zMod = 0; end
            
            plot(B.el, B.Z(:,k,i) - zMod, ...
                colors{mod(i-1,length(colors))+1}, ...
                'linewidth',2)
        end
        grid on;
        
        legend(num2str(1e-3*B.f(fIdx)'))
        xlabel('Elevation (deg)')
        ylabel('Magnitude (dB)')
        title(sprintf('Horizontal Beam Pattern @ %.1f degrees',B.el(k)))
        
    otherwise
        error('Unknown PLOTMODE parameter.  Should be one of ''vol'', ''surf'', or ''cont''.')
end

% get current axis handle
h = gca;


if ~isempty(SRC)
    if iscell(SRC)
        warning('SRC was entered, but could not determine which index to use!')
    else
        plot(SRC.az,SRC.el,'.')
    end
end
