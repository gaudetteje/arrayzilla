clc
clear
close all

% compute for different ranges
rng = 10;%[1 5 10];
N = numel(rng);

% create signal with AWGN, ~N(0,1)
ts.data = randn(2500,N);
ts.fs = 250e3;


%% test function
ts2 = az_armaloss(ts,rng,true);

%% plot data before/after filtering
cRange = [-60,60];
for k=1:N
    figure(1);
    title(sprintf('%d m - input',k))
    spectrogram(ts.data(:,k),64,62,256,ts.fs,'yaxis')
    set(gca,'CLim',cRange)
    colorbar

    figure(2);
    title(sprintf('%d m - output',k))
    spectrogram(real(ts2.data(:,k)),64,62,256,ts2.fs,'yaxis')
    set(gca,'CLim',cRange)
    colorbar
    
    tilefigs(2,2);
    pause
end