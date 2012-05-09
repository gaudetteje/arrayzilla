function x = genAudioTrack(cdata,wavfile,varargin)
% GENAUDIOTRACK  compiles echolocation data into an audio track
%
% genAudioTrack(CALLDATA,FNAME) takes in an Nx1 data struct containing
% reference channel data and reconstructs a "noise free" audio file, FNAME,
% with the correct time alignment.
% x = genAudioTrack(...) also returns the time series column vector
%
% The WAV file can be embedded into a beampattern video to show progression
% in time

% default parameters
D = 5;                         % slow audio playback by factor of D
tBuf = 0.1;                     % add short buffer to beginning
nbits = 16;


%% initialize full time series
N = numel(cdata);
tRef = cdata(1).t0 - tBuf;
tLen = tBuf + (cdata(N).t1(1) - cdata(1).t0(1));      % total duration
sLen = ceil(tLen * cdata(1).fs);

fprintf('Reconstructing original time series of length %g seconds and %d samples\n',tLen,sLen)
x = zeros(sLen,1);      % start with double precision, then reduce to int8 or int16 at end


%% iterate over each call
N = numel(cdata);                % total number of calls
for i = 1:N
    nSamp = size(cdata(i).data,1);
    
    t0 = cdata(i).t0 - tRef;    % find start time relative to x
    s0 = round(t0 * cdata(1).fs);        % find starting sample number in x
    idx = (s0:s0+nSamp-1);
    if isempty(idx), continue, end  % pass over empty blocks (failed during processing)
    
    % append data block to wav file with a taper to eliminate glitches in the audio
    x(idx) = raisedcos(numel(idx)) .* cdata(i).data;
end

% adjust maximum level to +/- 2^15 (for int16)
x = x ./ max(abs(x)/.999);  %(2^(nbits-1)

% write WAV file
wavwrite(x,cdata(1).fs/D,nbits,wavfile);

% plot data for sanity check
plot(x)
grid on
xlabel('Time (sec)')
ylabel('Amplitude (Volts)')

