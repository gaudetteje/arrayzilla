function plotTimeSeries(fname1,fname2,events)
% PLOTTIMESERIES  plots time series data from entire recording
%
% plotTimeSeries(FNAME1,FNAME2,EVENTMAP) plots the raw time series data for
% each call contained in EVENTMAP

N = numel(events);            % total number of events
ch1 = [72 1 110];
ch2 = [49 1];

figure
hold on;

% iterate over each call
for k = 1:N
    e = events(k);
    
    % sample index
    idx1 = (e.s0(1) : e.s1(1));
    idx2 = (e.s0(2) : e.s1(2));
    
    % read all data from selected channels
    data1 = read_SRZ(fname1,idx1,ch1) * (5 * 2^-16);
    for j=1:size(data1,2)
        data1(:,j) = data1(:,j) - mean(data1(:,j));
    end
    
    data2 = read_SRZ(fname2,idx2,ch2) * (5 * 2^-16);
    for j=1:size(data2,2)
        data2(:,j) = data2(:,j) - mean(data2(:,j));
    end
    
    % append to time index
    t1 = linspace(e.t0(1), e.t1(1), (e.s1(1)-e.s0(1))+1);
    t2 = linspace(e.t0(2), e.t1(2), (e.s1(2)-e.s0(2))+1);
    
    % plot current call
    plot(t1, data1)
    plot(t2, data2, '--')
    legend(num2str([ch1 ch2]'))

    clear idx* data* t*
end

grid on
xlabel('Time (sec)')
ylabel('Amplitude (Volts)')
