function res = az_filter(ts, src, a, varargin)
% AZ_FILTER  Realign multi-channel time series given the range to each sensor
%
% TS2 = az_filter(TS, SRC, A) takes in time series structure, TS, source
% location structure, SRC, and array structure, A, and returns a time
% series structure, TS2, aligned and optionally filtered around harmonics
% using FrFT-based instantaneous frequency estimation.


% optional inputs
DEBUG = true;  
FILTMODE = true;
REFMODE = 'fft';     % 'man'; %
nPad = 128;     % number of samples to keep before/after detected signal
nRef = 5;       % number of channels to use as reference
c = calcSoundSpeed(22.5);
if nargin > 2
    FILTMODE = varargin{1};
end

fprintf('\n\n***********************************************\n')


%% find best reference channel (or average IF of several channels together)
switch lower(REFMODE)
    case 'fft'
        % select the channel with the largest difference (yes, it's a crude SNR estimate)
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

        % find the strongest signals
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


%% compute delay based on ToA from source to array

% calculate time of arrival for each sensor
toa = src.rng / c;

% calculate time difference of arrival between each sensor
tdoa = toa-min(toa);

% calculate sample delay to closest whole sample
delta = round(tdoa .* ts.fs);


%% use reference channels to find start/stop time of signal

% filter out start/stop band (40-70kHz)
b = firls(80,[0 38e3 42e3 68e3 72e3 ts.fs/2]./(ts.fs/2),[0 0 1 1 0 0]);    %freqz(b0,1,4096,ts.fs)
Rhat = filtfilt(b,1,R);
Rhat = Rhat.^2;     % use signal energy

% find start/stop of each reference channel
for n=1:numel(refidx)
    t0(n) = find(Rhat(:,n) > 0.1*max(Rhat(:,n)),1,'first');
    t1(n) = find(Rhat(:,n) > 0.1*max(Rhat(:,n)),1,'last');
end
t0 = t0-nPad;    % to remove STFT transients
t1 = t1+nPad;    % to remove STFT transients

% correct for TDOA and average start/stop markers
t0 = floor(mean(t0 - delta(refidx)));
t1 = floor(mean(t1 - delta(refidx)));

%% realign and reduce data window

% initialize resulting time series data
nCh = size(ts.data,2);
nSamp = (t1-t0)+1;
%T = nSamp./ts.fs;       % calc pulse length
res.fs = ts.fs;
res.data = zeros(nSamp,nCh);

% remove extraneous samples before and after (with padding)
% iterate over each channel and realign relative to propagation delay
for ch = 1:nCh
    idx = (t0+delta(ch) : t1+delta(ch));
    res.data(:,ch) = ts.data(idx,ch) - mean(ts.data(idx,ch));
end


% option to skip filtering
if ~FILTMODE
    return
end


%% estimate IA and IF curve for reference channel(s)
H = hilbert(res.data);
for n = 1:numel(refidx)
    % estimate first component
    [IF1(:,n),p1(:,n)] = mca_ifestimate(H(:,n),1,.3,.3,.4);
    H1(:,n) = mca_iffilt(H(:,n),IF1(:,n));
    
    % estimate second component
    [IF2(:,n),p2(:,n)] = mca_ifestimate(H(:,n)-H1(:,n),1,.3,.3,.4);
    H2(:,n) = mca_iffilt(H(:,n),IF2(:,n));
    
end

% back out to expected starting time and take the average polynomial curve


% average IF polynomials
%IF{n} = IF{n} .* ts.fs;     % need to fix bug in mca_ifestimate when fs > 1




%% filter around IF estimate and separate harmonics on all channels




%debug
if DEBUG
    Fs = res.fs*1e-3;
    
    figure
    spectrogram(res.data(:, refch), 128, 120, 128, Fs, 'yaxis')
    xlabel('Time (ms)')
    ylabel('Frequency (kHz)')
    title(sprintf('Reference channel #%d',refch))
    hold on
    plot3((0:1/Fs:(length(IF)-1)/Fs), IF*1e-3,-50*ones(length(IF)),'k','linewidth',5)
    plot3((0:1/Fs:(length(IF)-1)/Fs), 2*IF*1e-3,-50*ones(length(IF)),'k','linewidth',5)
    plot3((0:1/Fs:(length(IF)-1)/Fs), 3*IF*1e-3,-50*ones(length(IF)),'k','linewidth',5)
    drawnow
end


% iterate over each channel to perform harmonic separation
fprintf('Filtering harmonic components\n');
for ch = 1:nCh

    % filter out 1st harmonic
    FM1 = mca_iffilt(res.data(:,ch), IF, ts.fs);

    % filter out 2nd harmonic
    FM2 = mca_iffilt(res.data(:,ch), 2.*IF, ts.fs);

    % filter out 3rd harmonic
    FM3 = mca_iffilt(res.data(:,ch), 3.*IF, ts.fs);

    % recombine time series for all channels
    res.data(:,ch) = real(FM1 + FM2 + FM3);
    
end

