%function [events, hdr] = az_detect(fname1,fname2)
% AZ_DETECT takes the raw binary recorder files and detects events in each
% channel.  Also removes bad channels from record/array definition
%
% events = az_detect(fname1,fname2)  returns a 1xM struct holding
% the start/stop sample numbers and times for M triggered events
%
% [events, hdr] = az_detect(fname1,fname2) also returns the header field
% information for both files


fprintf('\n\n***********************************************\n')
fprintf('Detecting triggered events in data files\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look for invalid CRC checksums (this is very time and memory consuming!)
%[stat, valid] = check_SRZ_CRC32(fname1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read header content 
fprintf('Reading header information from side 1...  \n\n')
%hdr(1) = read_header(fname1);

fprintf('\nReading header information from side 2...  \n\n')
%hdr(2) = read_header(fname2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify/correct data set alignment

% find start of continuous data blocks
hdr(1).block = [1; find(diff(hdr(1).count)~=1)+1];
hdr(2).block = [1; find(diff(hdr(2).count)~=1)+1];


% find indices to the first segment of each data block
N1 = numel(hdr(1).block);       % # blocks in side 1
fprintf('Found %d data blocks in "%s"\n',N1,fname1)
idx1 = zeros(N1+1,1);
for i=1:N1
    res = find(hdr(1).segment == hdr(1).block(i),1);
    assert(~isempty(res), 'Start of data block not aligned with any data segment!');      % THIS SHOULD NEVER HAPPEN
    idx1(i) = res;
end
idx1(end) = numel(hdr(1).segment)+1;  % append last segment index

N2 = numel(hdr(2).block);       % # blocks in side 2
fprintf('Found %d data blocks in "%s"\n',N2,fname2)
idx2 = zeros(N2+1,1);
for j=1:N2
    res = find(hdr(2).segment == hdr(2).block(j),1);
    assert(~isempty(res), 'Start of data block not aligned with any data segment!');      % THIS SHOULD NEVER HAPPEN
    idx2(j) = res;
end
idx2(end) = numel(hdr(2).segment)+1;  % append last segment index


%%%%%%%%%%%%%%%%%%%%
% first try looking for SMPTE data on auxiliary channel
%hdr(1).trig = hdr(1).trig - 2.5;
%hdr(2).trig = hdr(2).trig - 2.5;
%for i=1:N1
%    hdr(1).smpte{i} = SMPTE_dec(hdr(1).trig(idx1(i):idx1(i+1)-1), hdr(1).fs(1), 30);
%end

%%%%%%%%%%%%%%%%%%%%
% otherwise, cross-correlate data segments to match each block

% iterate over each data block to find matching pairs of segments
fprintf('\nSynchronizing data sequences...  ')
for i=1:N1
    blk1 = hdr(1).segment(idx1(i) : idx1(i+1)-1);
    for j=1:N2
        blk2 = hdr(2).segment(idx2(j) : idx2(j+1)-1);
        
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

[XC,lags] = xcorr(hdr(1).time(hdr(1).segment+1), hdr(2).time(hdr(2).segment+1));
[~,idx] = max(XC);
offset = lags(idx);
if offset
    warning('Offset between segments is not zero!  Possible synchronization issues between files.')
    if 1
        
        % correct segment offset
    end
end
fprintf('Done!\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = numel(hdr(1).segment);
events(N).s0 = 0;    % preallocate mem
events(N).s1 = 0;
events(N).t0 = 0;
events(N).t1 = 0;
for k = 1:N-1
    events(k).s0 = hdr(1).segment(k)+1;
    events(k).s1 = hdr(1).segment(k+1);
    events(k).t0 = hdr(1).time(events(k).s0);
    events(k).t1 = hdr(1).time(events(k).s1);
end
events(N).s0 = hdr(1).segment(N);
events(N).s1 = numel(hdr(1).count);
events(N).t0 = hdr(1).time(events(N).s0);
events(N).t1 = hdr(1).time(events(N).s1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Subfunction to read header info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function hdr = read_header(fname)
    res = read_SRZ_header_field(fname,[],[2 3 5 7 9]);
    
    % assign to hdr struct with appropriate names
    hdr.count = res(:,1);
    hdr.clock = res(:,2);
    hdr.gain = res(:,3);
    hdr.fs = res(:,4);
    hdr.trig = res(:,5);

    % unwrap clock cycles if overflow occurred
    k = find(diff(hdr.clock) < 1); k = [k; length(hdr.clock)];
    for kNum = 1:length(k)-1
        idx = k(kNum)+1:k(kNum+1);
        hdr.clock(idx) = hdr.clock(idx) + kNum*2^32;
    end
    
    % convert clock cycles to time
    hdr.time = (hdr.clock - hdr.clock(1)) ./ 106.25e6;
    
    % look for data segments by difference in clock time
    hdr.segment = [1; find(diff(hdr.clock) > 450)+1];
    fprintf('\n  Found %d data segments in "%s"\n', length(hdr.segment), fname);
