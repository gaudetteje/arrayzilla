function plotBeamPattern(B,varargin)
% PLOTBEAMPATTERN  plots an individual measured 3D beam
%
% Beam is a struct containing the 3-dimensional array of complex beam data,
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

fprintf('\n\n*****************************************\n')
fprintf('Plotting beam pattern\n\n')

% optional parameters
PLOTMODE = 'surf'; %'cont'
DATAMODE = 'fft'; %'hsa1'; %
switch (nargin)
    case 1
    case 2
        PLOTMODE = varargin{1};
    case 3
        PLOTMODE = varargin{1};
        DATAMODE = varargin{2};
    otherwise
        error('Incorrect number of input arguments')
end

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


%% hard coded parameters
smMETH = 'box';
smSIZE = 7;
B.Z = smoothn(B.Z,1e-4);                % interpolate NaN values (especially at missing corners)
B.Z = smooth3(B.Z,smMETH,smSIZE);       % smooth data in all 3 dimensions


% volume plotting options



% surface plotting options
cMap = 'hot';%flipud(hot);     %jet;    % colormap
dBrange = 35;           % colorscale depth
dBnorm = true;          % normalize to peak?
azView = 8;%20;            % azimuth angle
elView = 76;%20;            % elevation angle

% contour plotting
contourLev = -3;        % dB contour level for each frequency line [dB]
colors = {'k','b','g','m','c','r'};

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
        for i = 1:length(B.f)
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
        for i = 35 %1:length(B.f)
            
            % convert (az,el,rho) data points to cartesian coordinates
            z = -B.Z(:,:,i) .* cos(B.AZ*pi/180) .* cos(B.EL*pi/180);
            x = B.Z(:,:,i) .* sin(B.AZ*pi/180) .* cos(B.EL*pi/180);
            y = B.Z(:,:,i) .* sin(B.EL*pi/180);
            
            % find peak value
            dBpeak = max(max(B.Z(:,:,i)));
            fprintf('Peak @ %g kHz = %2.1f dB\n', B.f(i)*1e-3, dBpeak)
            
            % plot interpolated surface
            fh(i) = figure;%('MenuBar','none','ToolBar','none','NumberTitle','off');
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
            xlabel('az (\circ)','fontsize',16)
            ylabel('el (\circ)','fontsize',16)
            zlabel('dB','fontsize',16)
            
            set(gca,'fontsize',16);
            set(gcf,'color','w');
            
        end
        
        
        
    % plot 3D surface plots for each frequency bin
    case 'surf'
        
        % iterate over each frequency beam
        fh = nan(1,length(B.f));
        for i = 35 %1:length(B.f)
            
            % find peak value
            dBpeak = max(max(B.Z(:,:,i)));
            fprintf('Peak @ %g kHz = %2.1f dB\n', B.f(i)*1e-3, dBpeak)
            
            % plot interpolated surface
            fh(i) = figure;%('MenuBar','none','ToolBar','none','NumberTitle','off');
            surfc(B.AZ, B.EL, B.Z(:,:,i)-dBpeak);
            view(azView,elView)
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
            xlabel('az (\circ)','fontsize',16)
            ylabel('el (\circ)','fontsize',16)
            zlabel('dB','fontsize',16)
            
            set(gca,'fontsize',16);
            set(gcf,'color','w');
            
            
            pause(0.2)
        end
        %tilefigs(fh,2,3)
        
    % plot contours over all frequencies
    case 'cont'
        figure;
        hold on;
        dBlevels = contourLev * [1 1];
        for i = 1:length(B.f)
            dBpeak = max(max(B.Z(:,:,i)));
            fprintf('Peak @ %g kHz = %2.1f dB\n', B.f(i)*1e-3, dBpeak)
            contour(B.AZ, B.EL, B.Z(:,:,i),dBpeak+dBlevels,colors{mod(i-1,length(colors))+1},'linewidth',2);
        end
        grid on;
        legend(num2str(1e-3*B.f(:)))
        
    % plot horizontal slice
    case 'horz'
        k = 28; % select an elevation
        
        figure; hold on;
        for i = 1:length(B.f)
            % find peak value
            dBpeak = max(max(B.Z(:,:,i)));
            fprintf('Peak @ %g kHz = %2.1f dB\n', B.f(i)*1e-3, dBpeak)
            
            if dBnorm, zMod = dBpeak; else, zMod = 0; end
            
            plot(B.az, B.Z(k,:,i) - zMod, ...
                colors{mod(i-1,length(colors))+1}, ...
                'linewidth',2)
        end
        grid on;
        
        legend(num2str(1e-3*B.f(:)))
        xlabel('Azimuth (deg)')
        ylabel('Magnitude (dB)')
        title(sprintf('Horizontal Beam Pattern @ %.1f degrees',B.el(k)))
        
    % plot vertical slice
    case 'vert'
        k = 16; % select an azimuth
        
        figure; hold on;
        for i = 1:length(B.f)
            % find peak value
            dBpeak = max(max(B.Z(:,:,i)));
            fprintf('Peak @ %g kHz = %2.1f dB\n', B.f(i)*1e-3, dBpeak)
            
            if dBnorm, zMod = dBpeak; else, zMod = 0; end
            
            plot(B.el, B.Z(:,k,i) - zMod, ...
                colors{mod(i-1,length(colors))+1}, ...
                'linewidth',2)
        end
        grid on;
        
        legend(num2str(1e-3*B.f(:)))
        xlabel('Elevation (deg)')
        ylabel('Magnitude (dB)')
        title(sprintf('Horizontal Beam Pattern @ %.1f degrees',B.el(k)))
        
        
        
    otherwise
        error('Unknown PLOTMODE parameter.  Should be one of ''vol'', ''surf'', or ''cont''.')
end
