function [event, hdr] = az_detect_events(varargin)
% AZ_DETECT takes the raw binary recorder files and detects discontinuous
% events due to either trigger events, multiple recordings, or data corruption
%
% event = az_detect_events(FNAME1,FNAME2) returns a 1xM struct holding
%     the start/stop sample numbers and times for M triggered events
%
% event = az_detect_events(PREFIX) uses the string PREFIX to specify the
%     SRZ file pair
%
% [event, hdr] = az_detect_events(..) also returns the header field
%     information for both files
%
% event = az_detect_events(FNAME1,FNAME2,'eventfile') saves event struct
%     to the file 'eventfile.mat'
%
% event = az_detect_events(FNAME1,FNAME2,'eventfile','hdrfile') also
%     saves the header field information for all samples in 'hdrfile.mat'
%

fprintf('\n\n***********************************************\n')
fprintf('Detecting triggered events in data files\n')

switch nargin
    case 1
        prefix = varargin{1};
        if ~strcmp(prefix(end),'_')
            prefix(end+1) = '_';
        end
        fname1 = [prefix 'side1.srz'];
        fname2 = [prefix 'side2.srz'];
    case 2
        fname1 = varargin{1};
        fname2 = varargin{2};
    case 3
        fname1 = varargin{1};
        fname2 = varargin{2};
        eventfile = varargin{3};
    case 4
        fname1 = varargin{1};
        fname2 = varargin{2};
        eventfile = varargin{3};
        hdrfile = varargin{4};
    otherwise
        error('Incorrect number of input parameters')
end

% get filename prefix from current filename
if ~exist('prefix','var')
    prefix = regexp(fname1,'[_\-\ ]');
    prefix = fname1(1:prefix(end));
end

% assign eventfile name if not specified
if ~exist('eventfile','var')
    eventfile = [prefix 'event.mat'];
end

% assign hdrfile name if not specified
if ~exist('hdrfile','var')
    hdrfile = [prefix 'hdr.mat'];
end

% verify files exist
if ~exist(fname1,'file')
    error('AZ_DETECT_EVENTS:fnf', 'Could not locate file "%s"', fname1)
end
if ~exist(fname2,'file')
    error('AZ_DETECT_EVENTS:fnf', 'Could not locate file "%s"', fname2)
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
event(1:max(N)) = struct('s0',[0 0],'s1',[0 0],'t0',[0 0],'t1',[0 0]);
for n=1:max(N)
    event(n).eNum = n;     % assign an event number for reference

    % repeat for each side
    for m=1:2
        if (n < numel(hdr(m).event))
            event(n).s0(m) = hdr(m).event(idx(m,n));
            event(n).s1(m) = hdr(m).event(idx(m,n)+1)-1;
            event(n).t0(m) = hdr(m).time(event(n).s0(m));
            event(n).t1(m) = hdr(m).time(event(n).s1(m));
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save results to file
if ~isempty(eventfile)
    fprintf('\nSaving event index to "%s"...\n',eventfile)
    save(eventfile,'event');
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
    hdr.aux = res(:,5);     % Mx1 double
    
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
