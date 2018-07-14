
fname = '20130115_Dori2.AVI';

%% extract audio data from AVI file
[~,aud] = mmread(fullfile(pwd,fname),0);


%%
fs = aud.rate;
t = (0:size(aud.data,1)-1)./fs;
figure; plot(t,aud.data(:,1))

%%
idx = (6*fs:7*fs);
winlen = 128;
nfft = 512;
figure; spectrogram(aud.data(idx,1),hamming(winlen),winlen-2,nfft,fs,'yaxis')

