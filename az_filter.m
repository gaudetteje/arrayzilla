function res = az_filter(ts, src, a, varargin)
% AZ_FILTER  Realign multi-channel time series given the range to each sensor
%

% optional inputs
DEBUG = true;  
MODE = true;
if nargin > 2
    MODE = varargin{1};
end

fprintf('\n\n***********************************************\n')

DEBUG = true;

nPad = 128;     % number of samples to keep before/after detected signal

c = calcSoundSpeed(22);

% calculate time of arrival for each sensor
toa = src.rng / c;

% calculate time difference of arrival between each sensor
tdoa = toa-min(toa);

% calculate sample delay to closest whole sample
delta = round(tdoa .* ts.fs);

% find a good reference channel (ideally this would be based on SNR)
%fd = calc_spectrum(ts);        % use strongest signal in first harmonic band
%Fidx = (find(fd.freq > 20e3,1) : find(fd.freq < 50e3),-1);
%refch = max(sum(fd.magdb(:,Fidx),1))
%refch = find(tdoa == min(tdoa));
refch = 75;
refbd = 2;
refidx = az_chanindex(refch,refbd,a);
    

%% estimate start/stop time of strongest signal
% filter out start band (40-60kHz)
b0 = firls(80,[0 38e3 42e3 58e3 62e3 ts.fs/2]./(ts.fs/2),[0 0 1 1 0 0]);    %freqz(b0,1,4096,ts.fs)
ref0 = filtfilt(b0,1,ts.data(:,refidx));
ref0 = ref0.^2;     % use signal energy
t0 = find(ref0 > 0.1*max(ref0),1,'first');
t0 = t0-nPad;    % to remove STFT transients

% % filter out stop band (20-35kHz)
% b1 = firls(80,[0 18e3 22e3 33e3 37e3 ts.fs/2]./(ts.fs/2),[0 0 1 1 0 0]);    %freqz(b1,1,4096,ts.fs)
% ref1 = filtfilt(b1,1,ts.data(:,refidx));
% ref1 = ref1.^2;     % use signal energy
t1 = find(ref0 > 0.1*max(ref0),1,'last');
t1 = t1+nPad;    % to remove STFT transients

% initialize resulting time series data
nCh = size(ts.data,2);
nSamp = (t1-t0)+1;
%T = nSamp./ts.fs;       % calc pulse length
res.fs = ts.fs;
res.data = zeros(nSamp,nCh);


% iterate over each channel and realign relative to propagation delay
for ch = 1:nCh
    idx = (t0+delta(ch) : t1+delta(ch));
    res.data(:,ch) = ts.data(idx,ch) - mean(ts.data(idx,ch));
end

% option to skip filtering
if ~MODE
    return
end

% estimate IF of first harmonic using reference channel
IF = mca_ifestimate(hilbert(res.data(:, refidx))) .* res.fs;
%IF = mca_ifestimate(hilbert(res.data(nPad+1:end-nPad, refch)));%,res.fs);
%IF = [ones(nPad,1)*IF(1); IF; ones(nPad,1)*IF(end)];        % extend limits to padded buffers

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

% iterate over each channel
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

