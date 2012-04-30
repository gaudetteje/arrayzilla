function [events, hdr] = az_detect(fname1,fname2)
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
% read header content 
fprintf('Reading header information from side 1...  \n\n')
hdr(1) = read_header(fname1);

fprintf('\nReading header information from side 2...  \n\n')
hdr(2) = read_header(fname2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify/correct data set alignment, if sync data is available
if isfield(hdr,'sync')
    fprintf('\nSynchronizing data sequences...  \n')
    
    %  map synchronous trigger events to event index in both files
    [calls,idx(1,:),idx(2,:)] = intersect(hdr(1).sync, hdr(2).sync);
    N = numel(calls);
    M = numel(union(hdr(1).syn, hdr(2).sync));
    fprintf('  Found %d of %d synchronous trigger events.\n', N, M);
    
else
    % otherwise, just use first valid block for now and manually align data segments :(
    fprintf('\nNo synchronization data was found on auxiliary channels!  Assuming perfect alignment exists (check this manually).\n')
    
    % take smallest of 2 first blocks
    N = min([hdr(1).blockevent(2) hdr(2).blockevent(2)])-1;
    
    % assign fake call list and indices
    calls = 1:N;
    idx(1,:) = 1:N;
    idx(2,:) = 1:N;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assign event samples and times to struct
events(1:N) = struct('s0',[0 0],'s1',[0 0],'t0',[0 0],'t1',[0 0]);
for k=1:N
    events(k).s0(1) = hdr(1).event(idx(1,k));
    events(k).s1(1) = hdr(1).event(idx(1,k)+1)-1;
    events(k).s0(2) = hdr(2).event(idx(2,k));
    events(k).s1(2) = hdr(2).event(idx(2,k)+1)-1;

    events(k).t0(1) = hdr(1).time(events(k).s0(1));
    events(k).t1(1) = hdr(1).time(events(k).s1(1));
    events(k).t0(2) = hdr(2).time(events(k).s0(2));
    events(k).t1(2) = hdr(2).time(events(k).s1(2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Subfunction to read and process header info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = read_header(fname)
    res = read_SRZ_header_field(fname,[],[2 3 5 7 9]);      % can this result be converted directly to (u)int32?
    
    % assign to hdr struct with appropriate names
    hdr.count = res(:,1);   % Mx1 uint32
    hdr.clock = res(:,2);   % Mx1 uint32
    hdr.gain = res(:,3);    % 1x1 uint8
    hdr.fs = res(:,4);      % 1x1 double
    hdr.trig = res(:,5);    % Mx1 double
    
    % unwrap clock cycles if overflow occurred
    k = find(diff(hdr.clock) < 1); k = [k; length(hdr.clock)];
    for kNum = 1:length(k)-1
        idx = k(kNum)+1:k(kNum+1);
        hdr.clock(idx) = hdr.clock(idx) + kNum*2^32;
    end
    
    % convert clock cycles to time relative to start
    hdr.time = (hdr.clock - hdr.clock(1)) ./ 106.25e6;
    
    % look for discontinuous data blocks
    hdr.block = [1; find(diff(hdr.count)~=1)+1];
    fprintf('\n Found %d data blocks in "%s"\n',numel(hdr.block),fname)
    
    % look for trigger events by difference in clock time
    hdr.event = [1; find(diff(hdr.clock) > 450)+1];
    fprintf('\n  Found %d data segments (trigger events) in "%s"\n', length(hdr.event), fname);
    
    % find indices to the first segment of each data block
    hdr.blockevent = zeros(numel(hdr.block),1);
    for i=1:numel(hdr.block)
        res = find(hdr.event == hdr.block(i),1);
        assert(~isempty(res), 'Start of data block not aligned with any data segment!');      % THIS SHOULD NEVER HAPPEN
        hdr.blockevent(i) = res;
    end
    hdr.blockevent(end+1,1) = numel(hdr.event)+1;  % append last segment index

    % extract digital data from auxiliary channel, if possible
%     hdr.trig = hdr.trig - 2.5;                 % subtract DC offset
%     hdr.sync = zeros(size(hdr.event));
%     
%     % iterate over each event
%     for i=1:numel(hdr.event)
%         % compare level with 0 to recover DC level square wave
%         % synchronize to start pulse - remove prior data
%         % sample digital level every X seconds
%         % convert the digital serial word to a trigger count
%         hdr.sync(i) = res;
%     end
