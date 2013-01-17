function plotSpecArray(array,ts,varargin)
% PLOTSPECARRAY  displays a group of spectrograms from signals across the
% array face for debugging purposes and sanity checks

% plots spectrogram of a 5x5 set of channels around the array
chNum = [1 2 26 27 82 63 64 75 48 49 42 43 72 97 98 47 48 77 36 37 22 81 82 3 4];
bdNum = [1 1  2  2  2  1  1  2  2  2  1  1  1  2  2  1  1  1  2  2  1  1  1 2 2];

idx = az_chanindex(chNum, bdNum, array);

win = hamming(128);
nfft = 512;

fh = nan(1,length(idx));
for n=1:length(idx)
    fh(n) = figure('MenuBar','none','ToolBar','none','color','w');
    if ~idx(n), continue, end
    
    % zero pad data vector
    data = ts.data(:,idx(n));%[zeros(nfft,1); ts.data(:,idx(n)); zeros(nfft,1)];
    
    [~,F,T,P] = spectrogram(data, win,numel(win)-2,nfft,ts.fs*1e-3,'yaxis');
    surf(T,F,10*log10(abs(P)));
    shading interp
    view(0,90)
    %set(gca,'YDir','normal')
    
    set(gca,'YLim',[0 (ts.fs*1e-3)/2])
    set(gca,'XLim',[0 size(ts.data,1)/(ts.fs*1e-3)])
    % set(gca,'XLim',[nfft size(ts.data,1)+nfft]/(ts.fs*1e-3))
    set(gca,'CLim',[-80 -30]);
    colorbar
    
    title(sprintf('Channel %d, Side %d',chNum(n),bdNum(n)))
    xlabel('Time (ms)')
    ylabel('Frequency (kHz)')
    zlabel('Magnitude (dB)')
end

tilefigs(fh,5,5)
pause
close all