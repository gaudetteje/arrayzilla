clc
clear
close all

data = load('dori_1.mat');
sig = data.dori_1;
clear data

%% define parameters
fs=192e3;

alpha = 0.002;      % amplitude threshold - minimum envelope level for a valid event
gamma = 100;        % sample threshold - minimum number of samples between two continuous pulses

startpad = 100;
stoppad = 200;

min_dur = 500;      % throw away pulses having less than this number of samples

nfft = 256;
wind = hamming(64);
nover = length(wind)-4;

% debug flags
PLOT1 = false;
PLOT2 = false;
PLOT3 = false;
PLOT4 = true;

%% detect pulses based on instantaneous amplitude and signal duration
x=hilbert(sig);
y=find(abs(x) > alpha); %will give all the values of hilbert envelope above thresh
z=find(diff(y) > gamma); %gives intervals in samples between points; set this threshold to the minimum IPI of a pulse. 
call_start=[y(1); y(z+1)] - startpad;
call_stop=[y(z); y(end)] + stoppad;

% throw away calls that are too short in duration
valid_idx = find(call_stop-call_start > min_dur);
call_start = call_start(valid_idx);
call_stop = call_stop(valid_idx);

% verify start padding is not negative
call_start = max(call_start, 0);

% verify lengt of arrays
assert(length(call_start) == length(call_stop),'Received different number of start/stop points.  Check index construction.')

nCall = length(call_start);

%% plot time series with start/stop markers
if PLOT1
    close all
    plot(sig);
    grid on 
    hold on
    plot(abs(x),'--k');
    plot(call_start, zeros(size(call_start)), 'og'); %plotting the first sample above thresh of each call
    plot(call_stop, zeros(size(call_stop)), 'or'); %plotting the last sample above thresh of each call
end


%% build structure of individual data pulses
[pulse(1:nCall).data] = deal(zeros(500,1)); % initialize structure for speed
[pulse(1:nCall).time] = deal(zeros(500,1));
for i=1:nCall;
    pulse(i).data=sig(call_start(i):call_stop(i));
    pulse(i).time=(call_start(i):call_stop(i))'./ fs;
end


%% plot spectrograms of data
if PLOT2
    for j=1:length(pulse);
        spectrogram(pulse(j).data, wind, nover, nfft, fs, 'yaxis')
        pause(0.3)
    end
end

%% iterate Hilbert spectral analysis (HSA) over each call
tic
[pulse(1:nCall).IA1] = deal(zeros(500,1));
[pulse(1:nCall).IF1] = deal(zeros(500,1));
[pulse(1:nCall).IA2] = deal(zeros(500,1));
[pulse(1:nCall).IF2] = deal(zeros(500,1));

% add 'true' as a last parameter in any of the following mca_x() functions to plot the processing steps - useful for validation of results
for n=1:nCall;
    
    % process first harmonic (CAUTION:  this is not necessarily FM1)
    pulse(n).if1 = mca_ifestimate(pulse(n).data,1,0.5,0.5,0.3,2);               %estimate of the instantaneous frequency of strongest harmonic: eq B5 in DiCecco 2013
    pulse(n).fm1 = mca_iffilt(pulse(n).data,pulse(n).if1,1);                              %demodulating, filtering and extracting component as in fig. 3.
    [~,pulse(n).IA1,pulse(n).IF1] = mca_extract(pulse(n).fm1,fs,0,0,1);           %accurate estimate of instant freq and amplitude
    
    % process second harmonic
    pulse(n).if2 = mca_ifestimate(pulse(n).data-real(pulse(n).fm1),1,0.5,0.5,0.3,2);  % subtract FM1 from original time series and repeat processing
    pulse(n).fm2 = mca_iffilt(pulse(n).data-pulse(n).fm1,pulse(n).if2,1);
    [~,pulse(n).IA2,pulse(n).IF2] = mca_extract(pulse(n).fm2,fs,0,0,1);

end
toc

%% plot spectrograms and overlaid HSA results for each call
if PLOT3
    close all
    figure(1)
    figure(2)
    figure(3)
    tilefigs(2,3)
    
    for n=1:nCall
        figure(1)%,'NextPlot','new')
        [S,F,T,P] = spectrogram(pulse(n).data, wind, nover, nfft, fs*1e-3, 'yaxis');
        ah = pcolor(T,F,10*log10(abs(P)));
        set(ah,'LineStyle','none')
        ax = axis;
        axis(ax);
        cLim = get(gca,'CLim');
        cLim(1) = cLim(2)-50;
        set(gca,'CLim',cLim);
        colorbar
        hold on
        plot((pulse(n).time - pulse(n).time(1))*1e3, pulse(n).if1*fs*1e-3, '--m')
        plot((pulse(n).time - pulse(n).time(1))*1e3, pulse(n).if2*fs*1e-3, '--m')
        hold off
        drawnow
        
        figure(2)
        [S,F,T,P] = spectrogram(real(pulse(n).fm1), wind, nover, nfft, fs*1e-3, 'yaxis');
        ah = pcolor(T,F,10*log10(abs(P)));
        set(ah,'LineStyle','none')
        axis(ax);
        set(gca,'CLim',cLim);
        colorbar
        hold on
        plot((pulse(n).time - pulse(n).time(1))*1e3, pulse(n).IF1*1e-3, '--m')
        hold off
        drawnow
        
        figure(3)
        [S,F,T,P] = spectrogram(real(pulse(n).fm2), wind, nover, nfft, fs*1e-3, 'yaxis');
        ah = pcolor(T,F,10*log10(abs(P)));
        set(ah,'LineStyle','none')
        axis(ax);
        set(gca,'CLim',cLim);
        colorbar
        hold on
        plot((pulse(n).time - pulse(n).time(1))*1e3, pulse(n).IF2*1e-3, '--m')
        hold off
        drawnow
        
        pause(1)
    end
end


%% plot IF/IA over time

if PLOT4
    cRange = [-60 -20];
    cMap = jet;
    for n=1:nCall
        % plot first harmonic (NOTE: ccplot requires row vectors)
        lh1 = ccplot(1e3*pulse(n).time', 1e-3*pulse(n).IF1',db(pulse(n).IA1)', cMap, cRange);
        hold on;
        
        % plot second harmonic
        lh2 = ccplot(1e3*pulse(n).time', 1e-3*pulse(n).IF2',db(pulse(n).IA2)', cMap, cRange);
    end
    colorbar
end
