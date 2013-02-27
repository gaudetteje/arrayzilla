function [events, hdr] = az_detect_events(fname1,fname2,varargin)
% AZ_DETECT takes the raw binary recorder files and detects discontinuous
% events due to either trigger events, multiple recordings, or data corruption
%
% events = az_detect_events(fname1,fname2)  returns a 1xM struct holding
%     the start/stop sample numbers and times for M triggered events
% [events, hdr] = az_detect_events(fname1,fname2) also returns the header
%     field information for both files
% events = az_detect_events(fname1,fname2,'eventfile')  saves events struct
%     to the file 'eventfile.mat'.
% events = az_detect_events(fname1,fname2,'eventfile','hdrfile')  also
%     saves the header field information for all samples in 'hdrfile.mat'.

fprintf('\n\n***********************************************\n')
fprintf('Detecting triggered events in data files\n')

if ~exist(fname1,'file')
    error('AZ_DETECT:fnf', 'Could not locate file "%s"', fname1)
end
if ~exist(fname2,'file')
    error('AZ_DETECT:fnf', 'Could not locate file "%s"', fname2)
end

hdrfile = [];
switch nargin
    case 2
        prefix = regexp(fname1,'[_\-\ ]');      % use current filename
        eventfile = [fname1(1:prefix(end)) 'events.mat'];
    case 3
        eventfile = varargin{1};
    case 4
        eventfile = varargin{1};
        hdrfile = varargin{2};
    otherwise
        error('Incorrect number of input parameters')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read header content 
fprintf('\nReading header information from side 1...  \n')
hdr(1) = read_header(fname1);

fprintf('\nReading header information from side 2...  \n')
hdr(2) = read_header(fname2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify/correct data set alignment, if sync data is available
if isfield(hdr,'sync')
    fprintf('\nSynchronizing data sequences...  \n')
    
    %  match synchronous trigger events to event index in both files
    [match,idx(1,:),idx(2,:)] = intersect(hdr(1).sync, hdr(2).sync);
    N = numel(match);
    M = numel(union(hdr(1).syn, hdr(2).sync));
    fprintf('  Found %d of %d synchronous trigger events.\n', N, M);
    
else
    % otherwise, just use first valid block for now and manually align data segments :(
    fprintf('\nNo synchronization data was found on auxiliary channels!  Assuming perfect alignment exists\n\t(correct this manually by adjusting event index).\n')
    
    % determine number of events in each side
    N = [numel(hdr(1).event) numel(hdr(2).event)]-1;
    
    % assign event list and indices
    idx = [1:N(1) nan(1,N(2)-N(1)) ; 1:N(2) nan(1,N(1)-N(2))];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assign event samples and times to struct
events(1:max(N)) = struct('s0',[0 0],'s1',[0 0],'t0',[0 0],'t1',[0 0]);
for k=1:max(N)
    events(k).eNum = k;     % assign an event number for reference
    
    if (k < numel(hdr(1).event))
        events(k).s0(1) = hdr(1).event(idx(1,k));
        events(k).s1(1) = hdr(1).event(idx(1,k)+1)-1;
        events(k).t0(1) = hdr(1).time(events(k).s0(1));
        events(k).t1(1) = hdr(1).time(events(k).s1(1));
    end
    if (k < numel(hdr(2).event))
        events(k).s0(2) = hdr(2).event(idx(2,k));
        events(k).s1(2) = hdr(2).event(idx(2,k)+1)-1;
        events(k).t0(2) = hdr(2).time(events(k).s0(2));
        events(k).t1(2) = hdr(2).time(events(k).s1(2));
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save results to file
if ~isempty(eventfile)
    fprintf('\nSaving event index to "%s"...\n',eventfile)
    save(eventfile,'events');
end
if ~isempty(hdrfile)
    fprintf('Saving header information to "%s"...\n',hdrfile)
    save(hdrfile,'hdr');
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
    hdr.aux = res(:,5);    % Mx1 double
    
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
    fprintf('Found %d data blocks in "%s"\n',numel(hdr.block),fname)
    
    % look for trigger events by difference in clock time - append last sample
    hdr.event = [1; find(diff(hdr.clock) > 450)+1; numel(hdr.count)+1];
    fprintf('Found %d data segments (trigger events) in "%s"\n', length(hdr.event)-1, fname);
    
    % find index to the first event in each data block
    hdr.blockidx = zeros(numel(hdr.block),1);
    for i=1:numel(hdr.block)
        res = find(hdr.event == hdr.block(i),1);
        assert(~isempty(res), 'Start of data block not aligned with any data segment!');      % THIS SHOULD NEVER HAPPEN
        hdr.blockidx(i) = res;
    end
    hdr.blockidx(end+1,1) = numel(hdr.event);  % append last segment index
    
    %%% use FFT to look for 10kHz digital data stream on auxiliary channel
    % extract digital data from auxiliary channel, if possible
%     hdr.aux = hdr.aux - mean(hdr.aux);                 % subtract DC offset
%     % evaluate sync for 
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
