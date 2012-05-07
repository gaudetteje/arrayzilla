function res = az_filter(ts, src, a, varargin)
% AZ_FILTER  Realign multi-channel time series given the range to each sensor
%
% TS2 = az_filter(TS, SRC, A) takes in time series structure, TS, source
% location structure, SRC, and array structure, A, and returns a time
% series structure, TS2, aligned and optionally filtered around harmonics
% using FrFT-based instantaneous frequency estimation.
%
% TS2 is a time series struct if data is aligned without harmonic
% separation.  Otherwise, TS2 is a cell of time series structs holding each
% separated harmonic.
%
% Sequence of automated processing events:
% - find a set of reference channels with high SNR
% - calculate time of arrival from localized source to each mic channel
% - estimate global start/stop time from bandlimited energy in reference channels
% - realign data matrix for all channels
% +++ if filtering mode enabled:
% - estimate average IF for each harmonic from reference channels
% - separate harmonics on all channels

% optional inputs
DEBUG = false;
FILTMODE = true;
REFMODE = 'fft';     % 'man'; %
nPad = 128;     % number of samples to keep before/after detected signal
nRef = 1;       % number of channels to use as reference
BW = .035;      % lowpass filter bandwidth (normalized) for mca_iffilt
c = calcSoundSpeed(22.5);
if nargin > 2
    FILTMODE = varargin{1};
end

fprintf('\n\n***********************************************\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find best reference channel(s) in data set
switch lower(REFMODE)
    case 'fft'
        % select the channel with the largest SNR difference
        fd = calc_spectrum(ts,size(ts.data,1),@rectwin,1);        % compute spectrum
        freq = fd.freq;
        mag = sgolayfilt(fd.mag,1,101);
        clear fd;

        % sum in-band energy between 30kHz and 80kHz
        fl0 = find(freq > 30e3,1);
        fh0 = find(freq < 80e3,1,'last');
        E0 = sum(mag(fl0:fh0,:),1);

        % sum out of band energy outside 25kHz and 85kHz
        fl1 = find(freq < 25e3,1);
        fh1 = find(freq < 85e3,1,'last');
        N0 = sum(mag([1:fl1 fh1:end],:),1);

        % find the strongest signals (yes, it's crude but it works!)
        [~,idx] = sort(E0./N0,'descend');
        refidx = idx(1:nRef);
        
    case 'man'
        % alternatively, we can manually select references
        refch = 75;
        refbd = 2;
        refidx = az_chanindex(refch,refbd,a);
        
    otherwise
        error('Unknown filtering mode')
end

R = ts.data(:,refidx);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute delay based on ToA from source to array

% calculate time of arrival for each sensor
toa = src.rng / c;

% calculate time difference of arrival between each sensor
tdoa = toa-min(toa);

% calculate sample delay to closest whole sample
delta = round(tdoa .* ts.fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use reference channels to find relative start/stop time of signal

% filter out start/stop band (40-70kHz)
b = firls(80,[0 38e3 42e3 68e3 72e3 ts.fs/2]./(ts.fs/2),[0 0 1 1 0 0]);    %freqz(b0,1,4096,ts.fs)
Rhat = filtfilt(b,1,R);
Rhat = Rhat.^2;     % use signal energy

% find start/stop of each reference channel
for n=1:numel(refidx)
    t0(n) = find(Rhat(:,n) > 0.1*max(Rhat(:,n)),1,'first');
    t1(n) = find(Rhat(:,n) > 0.1*max(Rhat(:,n)),1,'last');
end

% correct for TDOA and take the 95% (mean - 3 std. dev.) start/stop markers
t0 = max(floor(mean(t0-delta(refidx)) - 3*std(t0-delta(refidx))) - nPad, 1);
t1 = min(ceil(mean(t1-delta(refidx)) + 3*std(t1-delta(refidx))) + nPad, size(R,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% realign and reduce data window

% initialize resulting time series data
nCh = size(ts.data,2);
nSamp = (t1-t0)+1;
T = nSamp./ts.fs;       % calc pulse length
fprintf('Call duration is %.3g ms\n', 1e3*(T-(2*nPad/ts.fs)))

% iterate over each channel and realign relative to receive delay (with padding)
res.fs = ts.fs;
res.data = zeros(nSamp,nCh);
for ch = 1:nCh
    idx = (t0+delta(ch) : t1+delta(ch));        % MIGHT NEED ASSERTION TO PREVENT OUT OF BOUNDS ERROR
    res.data(:,ch) = ts.data(idx,ch) - mean(ts.data(idx,ch));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate IA and IF curve for reference channel(s)

% bailout option to skip filtering
if ~FILTMODE, return, end

% convert real data to its analytic form
H = hilbert(res.data);

% iterate over each reference channel
for n = 1:numel(refidx)
    % estimate first component
    IF1 = mca_ifestimate(H(:,refidx(n)),1,.4,.4,.4);
    H1 = mca_iffilt(H(:,refidx(n)),IF1,1);
    
    % estimate second component
    IF2 = mca_ifestimate(H(:,refidx(n))-H1,1,.4,.4,.4);
    
    % swap harmonics if necessary
    if sum(IF1) > sum(IF2)
        FM1(:,n) = IF2;
        FM2(:,n) = IF1;
    else
        FM1(:,n) = IF1;
        FM2(:,n) = IF2;
    end
    
    % plot spectrogram for debugging
    if DEBUG
        fh0 = figure; spectrogram(real(H(:,refidx(n))),64,62,256,res.fs,'yaxis');
        hold on;
        clim = get(gca,'clim');
        set(gca,'clim',[clim(2)-50 clim(2)]);
        
        plot3((1:length(FM1))./res.fs, FM1.*res.fs, clim(2)*ones(length(FM1)), 'g', 'linewidth',2)
        plot3((1:length(FM1))./res.fs, (FM1+BW).*res.fs, clim(2)*ones(length(FM1)), '--g', 'linewidth',1)
        plot3((1:length(FM1))./res.fs, (FM1-BW).*res.fs, clim(2)*ones(length(FM1)), '--g', 'linewidth',1)
        
        plot3((1:length(FM2))./res.fs, FM2.*res.fs, clim(2)*ones(length(FM2)), 'm', 'linewidth',2)
        plot3((1:length(FM2))./res.fs, (FM2+BW).*res.fs, clim(2)*ones(length(FM2)), '--m', 'linewidth',1)
        plot3((1:length(FM2))./res.fs, (FM2-BW).*res.fs, clim(2)*ones(length(FM2)), '--m', 'linewidth',1)
        
        title(sprintf('Idx = %d, Bd = %d, Ch = %d, (x = %gm, y = %gm)', refidx(n), a.bd(refidx(n)), a.ch(refidx(n)), a.xPos(refidx(n)), a.yPos(refidx(n))))
    end
end

% take the average IF curves
FM1 = mean(FM1,2);
FM2 = mean(FM2,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% filter around IF estimate and separate harmonics on all channels

% design flat cheby II filter
[b,a] = cheby2(8,40,.125);

% iterate over each channel to perform harmonic separation
fprintf('Filtering harmonic components\n');
%res.data = [];
res.fm1 = zeros(size(res.data));
res.fm2 = zeros(size(res.data));
for n = 1:nCh

    % filter out 1st harmonic
    res.fm1(:,n) = mca_iffilt(res.data(:,n), FM1, 1, b, a);

    % filter out 2nd harmonic
    res.fm2(:,n) = mca_iffilt(res.data(:,n), FM2, 1, b, a);

    % recombine into "filtered" data set
    res.data(:,n) = res.fm1(:,n) + res.fm2(:,n);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% debug
if 0 %DEBUG
    clim = [-110 -60];
    for n=1:10:209
        
        fh1 = figure; spectrogram(real(res.data(:,n)) ,64,62,256,ts.fs,'yaxis'); colorbar; set(gca,'clim',clim)
        fh2 = figure; spectrogram(real(res.fm1(:,n)),64,62,256,ts.fs,'yaxis'); colorbar; set(gca,'clim',clim)
        fh3 = figure; spectrogram(real(res.fm2(:,n)) ,64,62,256,ts.fs,'yaxis'); colorbar; set(gca,'clim',clim)
        fh4 = figure; spectrogram(real(res.data(:,n)) - real(res.fm1(:,n)) - real(res.fm2(:,n)),64,62,256,ts.fs,'yaxis'); colorbar; set(gca,'clim',clim)
        hold on
        
        tilefigs([fh1 fh2 fh3 fh4],2,2)
        drawnow
        pause
        
        %close([fh1 fh2 fh3 fh4])
    end
end
