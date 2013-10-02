
% create signal of AWGN, ~N(0,1)
ts.data = randn(2500,1);
ts.fs = 250e3;

ts2 = az_armaloss(ts,10);

%% plot data before/after filtering
figure;
spectrogram(ts.data,64,62,256,ts.fs,'yaxis')

