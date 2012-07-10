function calibrate_sensors(varargin)
% CALIBRATE_SENSORS  loads calibration data from microphones and iterates
% over each microphone to generate an ARMA filter response
%
% calibrate_sensors(fname1,fname2) loads data in each file and finds ARMA
% transfer function to correct for differences in amplitude/phase response
%
% calibrate_sensors(fname1,fname2,calfile,CH) updates the calibration data
% in calfile only for channels in the array CH.

FORCEDET = false;

fname1 = '/Users/jasongaudette/Documents/MATLAB/arrayzilla/calibration_data/20111104_calibration/caldata_side1.srz';
fname2 = '/Users/jasongaudette/Documents/MATLAB/arrayzilla/calibration_data/20111104_calibration/caldata_side2.srz';
calfile = '/Users/jasongaudette/Documents/MATLAB/arrayzilla/calibration_data/20111104_calibration/caldata_bk1.bin';

if ~exist('TDOA_frame','file')
    addpath('~/src/simmons_svn/flightroom_tools/primary_analysis/');
end
if ~exist('read_arraybin','file')
    addpath('~/src/simmons_svn/FieldArrayTracker/');
end

% % prompt for filename if not entered
% switch nargin
%     case {0,1}
%         
%         % prompt user for data file
%         [fname, pname] = uigetfile( { ...
%             '*.srz','Recorder Data (*.srz)'; ...
%             '*.*',  'All Files (*.*)'}, ...
%             'Select Side 1 data file', ...
%             'MultiSelect', 'on');
%         
%         if iscell(fname) && length(fname) == 2
%             % assign 2nd filename if multiple selected
%             fname1 = fullfile(pname,fname{1});
%             fname2 = fullfile(pname,fname{2});
%         else
%             % otherwise get 2nd file name
%             fname1 = fullfile(pname,fname);
%             [fname, pname] = uigetfile( { ...
%                 '*.srz','Recorder Data (*.srz)'; ...
%                 '*.*',  'All Files (*.*)'}, ...
%                 'Select Side 1 data file');
%             fname2 = fullfile(pname,fname);
%         end
%         clear fname pname
%         
%     case 2
%         % assign filenames
%         fname1 = varargin{1};
%         fname2 = varargin{2};
%         
%     otherwise
%         error('Bad number of input parameters.  Better luck next time!');
% end

%%%%%%%%%%%%%%%%

% identify known bad cal trials
%TBD


% load reference signal from B&K mic (relative to input @ microphone)
ref.data = read_arraybin(calfile,50)';
ref.avg = sum(ref.data,2);
ref.avg = ref.avg - mean(ref.avg);
ref.fs = 5e5;

ref.env = abs(hilbert(ref.avg));
ref.s0 = max(find(ref.env > .1*max(ref.env), 1, 'first') - 100, 0);
ref.s1 = find(ref.env > .1*max(ref.env), 1, 'last') + 100;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define physical array parameters
array = az_positions([19 12],[4 5]*.0254,[9/8 9/10]);

% assign channel/board mapping to array struct
array = az_channelmap(array);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detect calls in data files and return timestamps
prefix = regexp(fname1,'[_\-\ ]');      % use current filename
callmap = [fname1(1:prefix(end)) 'callmap.mat'];
if exist(callmap,'file') && ~FORCEDET
    fprintf('\nCall index file already exists!  Loading call data in "%s"...\n\n',callmap);
    load(callmap,'callmap');
else
    callmap = az_detect(fname1,fname2);
    fprintf('\nSaving call index to "%s"...\n\n',callmap)
    save(callmap,'callmap');
end


% select threshold for minimum channel gain
thresh = 1e4;
N = 700;  % length of signals
M = numel(array.bd);    % number of array elements

% optimal rotation for fractional fourier transform of pulse
a = .69;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iterate over each detected data segment
ch = zeros(1,2276);
val = zeros(1,2276);
for k=120:2276 %numel(callmap)
    tic
    % segment index
    idx = (callmap(k).s0 : callmap(k).s0+N-1);
    
    % load channel signal
    res1 = read_SRZ(fname1,idx,1:112);
    res2 = read_SRZ(fname2,idx,1:112);
    res = [res1 res2]-2^15;
    clear res1 res2

    % iterate over each channel to find signal with largest fractional peak
    z0 = zeros(1,M);
    for m = 1:M
        RWT = abs(frft(hilbert(res(:,m)),a));
        z0(m) = max(RWT);
    end
    
    [val(k),ch(k)] = max(z0);
    fprintf('\n[%d] CH = %d (%g)',k,ch(k),val(k));
    
    if val(k) < thresh
        fprintf('  Possible bad channel...')
    end
    
    toc
    %% calculate system impulse response
    
%     % Take Fourier transform of input and output signals
%     nfft = 8192;
%     X = fft(ts1.data,nfft);
%     Y = fft(ts2.data,nfft);
%     F = (0:(nfft-1))*(ts1.fs/nfft);
%     
%     % Compute the frequency response at all frequencies
%     H = Y./X;
%     
%     % Find which fft samples are in the chirp freq range
%     f1 = 10e3; % chirp low freq
%     f2 = 100e3; % chirp high freq
%     k1 = ceil((f1/ts1.fs)*nfft)+1;  % which fft sample corresponds to f1
%     k2 = floor((f2/ts1.fs)*nfft)+1; % which fft sample corresponds to f2
%     
%     % Set H to -30dB below f1 and -40dB above f2 (lowpass filter effect)
%     H(1:(k1-1)) = 10.^(-30/20).*ones(k1-1,1);
%     H((k2+1):(nfft/2)) = 10.^(-40/20).*ones((nfft/2)-k2,1);
%     H((nfft/2)+1) = 0.1;
%     H(((nfft/2)+2):nfft) = conj(flipud(H(2:(nfft/2))));
%     
%     % Smooth output signal's magnitude response, phase is effectively 0
%     % Not sure why there are lots of deviations:
%     %   Could be due to reverb or other reflections interfering with
%     %   recorded microphone signal)
%     Hdmag = abs(H);
%     Hdmag = sgolayfilt(Hdmag,3,101);       % use SG filter
%     Hdmag = filtfilt(ones(50,1)./50,1,Hdmag);        % use MA filter
%     [Hdr,Hdi] = pol2cart(zeros(size(Hdmag)),Hdmag);
%     Hd = Hdr + 1i.*Hdi;
%     
%     if PLOTFLAG
%         
%         % plot input output response
%         figure(199);
%         subplot(2,1,1)
%         %plot(F,db(abs(X)),F,db(abs(Y)),'r',F,db(abs(H)),'k',F,db(abs(Hd)),'--k')
%         plot(F,db(abs(H)),'DisplayName',fname);
%         grid on
%         hold on
%         %lh = legend(fname);  %'Input Signal','Output Signal')
%         %set(lh,'Interpreter','none');
%         ylabel('Magnitued [dB]')
%         title('Magnitude response of Input/Output')
% 
%         subplot(2,1,2)
%         plot(F,unwrap(angle(H))*pi/180,'DisplayName',fname)
%         %plot(F,unwrap(angle(X))*pi/180,F,unwrap(angle(Y))*pi/180,'r',F,unwrap(angle(H))*pi/180,'k')
%         grid on
%         hold on
%         %set(gca,'DisplayName',fname)
%         title('Phase response of Input/Output')
%         ylabel('Phase [Degrees]')
%         xlabel('Frequency [Hz]')
%         
%     end
end

% save updated microphone calibration data



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%     Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,z] = find_peaks(Z,x,y,varargin)
% FIND_PEAKS  find peak locations in 2D surface
%
% [xMax,yMax,zMax] = FIND_PEAKS(Z,x,y) returns the highest peak coordinates
%     and value given the surface, Z, and index vectors, x and y.
% [x,y,z] = FIND_PEAKS(Z,x,y,N) returns column vectors of the N highest peaks

N = 1;
if nargin > 3
    N = varargin{1};
end

% sort surface as column and save index
Zcol = Z(:);
[val,idx] = sort(Zcol,'descend');

% extract N highest peaks
[ypk,xpk] = ind2sub(size(Z),idx(1:N));
x = x(xpk);
y = y(ypk);
z = val(1:N);