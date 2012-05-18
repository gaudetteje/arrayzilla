%testbed for data sync

% use aligndata from Comms toolbox?

% iterate over each data block to find matching pairs of segments

% attempt to cross-correlate data segments to match each data block
fprintf('\nSynchronizing data sequences...  ')
for i=1:N1
    blk1 = hdr(1).event(idx1(i) : idx1(i+1)-1);
    for j=1:N2
        blk2 = hdr(2).event(idx2(j) : idx2(j+1)-1);
        
        % compare blocks until match is found
        % iterate over each possible delay (+/-)
        minlag = min(numel(blk1), numel(blk2));
        %maxlag = max(numel(blk1), numel(blk2));
        %lags = (-maxlag:maxlag);
        %[~,lagidx] = sort(lags);
        %lags = lags(lagidx);
        
        
        %for k = lags
            % look forward
        %    sum(blk1(1:minlag) - blk2(k:)
            
            % look backward
        %    sum(blk1(1:minlag))
        %end
        
        %[XC,lag] = xcorr(diff(blk1(1:minlag)),diff(blk2(1:minlag)));
        %[~,lagIdx] = max(XC);
        %lag = lag(lagIdx);
        
        if lag > 0
            
        elseif lag < 0 
            
        else
            
        end
        
        
        plot(XC)
        title(sprintf('Lag = %d',lag))
        
        
        % use MF on diff?
        
        if 0
            % remove overlapping segments from blk1
            break
        end
    end
    % assign matching pairs if found
end

[XC,lags] = xcorr(hdr(1).time(hdr(1).event+1), hdr(2).time(hdr(2).event+1));
[~,idx] = max(XC);
offset = lags(idx);
if offset
    warning('Offset between segments is not zero!  Possible synchronization issues between files.')
    if 1
        
        % correct segment offset
    end
end
