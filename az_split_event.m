function [call, ts] = az_split_event(varargin)
% AZ_SPLIT_EVENT searches raw binary data for one or more calls in each event
%
% CALL = az_split_event(FNAME1,FNAME2,EVENT) returns a CALL index struct
%       from data in the SRZ file pair (FNAME1, FNAME2) referenced in
%       EVENT struct
%
% CALL = az_split_event(PREFIX,EVENT) uses the string PREFIX to specify the
%       SRZ file pair
%
% CALL = az_split_event(FNAME1,FNAME2,EVENT,ARRAY) uses the channel
%       mapping in ARRAY to eliminate dead channels from analysis
%
% CALL = az_split_event(...,CALLFILE) entering a string as the last
%       parameter will save the resulting call struct to CALLFILE
%
% [CALL,ARRAY] = az_split_event(...) also returns the ARRAY struct with a
%       list of detected bad channels
%
% The function returns a new Mx1 structure, CALLS, with identical fields
%     as EVENT.  Long events (100ms or greater) are split into multiple
%     overlapping events before searching for calls.  This is to avoid
%     memory issues and allows processing continuous recorded data where no
%     trigger events are found.
%
% Note:  EVENT is generated using az_detect_events and ARRAY is generated
%     using az_positions and az_channelmap.

% TODO:
% - detect bad channels using sum(E,1) to look across time; return list

fprintf('\n***********************************************\n')

% set default parameters
BLOCKSIZE = 23666;              % process blocks of 100ms maximum

gamma = 1e4;                   % normalized threshold for amplitude detection
%nCh = 5;                        % number of channels required for threshold
nPad = 200;                     % number of samples to pad around detected calls
DEBUG = false;

% init parameters
array = [];
cNum = 0;                       % counter for number of detected calls

% handle optional inputs
switch nargin
    case 2
        prefix = varargin{1};
        event = varargin{2};
        if ~strcmp(prefix(end),'_')
            prefix(end+1) = '_';
        end
        fname1 = [prefix 'side1.srz'];
        fname2 = [prefix 'side2.srz'];
    case 3
        fname1 = varargin{1};
        fname2 = varargin{2};
        event = varargin{3};
    case 4
        fname1 = varargin{1};
        fname2 = varargin{2};
        event = varargin{3};
        res = varargin{4};
        if ischar(res)
            callfile = res;
        else
            array = res;
        end
    case 5
        fname1 = varargin{1};
        fname2 = varargin{2};
        event = varargin{3};
        array = varargin{4};
        callfile = varargin{5};
    otherwise
        error('Incorrect number of parameters entered')
end

% get filename prefix from current filename
if ~exist('prefix','var')
    prefix = regexp(fname1,'[_\-\ ]');
    prefix = fname1(1:prefix(end));
end

% assign callfile name if not specified
if ~exist('callfile','var')
    callfile = [prefix 'call.mat'];
end

% verify files exist
if ~exist(fname1,'file')
    error('AZ_SPLIT_EVENT:fnf', 'Could not locate file "%s"', fname1)
end
if ~exist(fname2,'file')
    error('AZ_SPLIT_EVENT:fnf', 'Could not locate file "%s"', fname2)
end


%% load a valid event struct from file if a filename was entered
if ischar(event)
    try
        load(event,'event')
    catch
        error('Event file not found:  "%s"',event)
    end
end

% verify event is a structure and has appropriate fields
if ~exist('event','var') || ~isstruct(event)
    error('Event file "%s" does not contain "event" struct!',event)
end
fields = fieldnames(event);
if ~all(ismember({'s0','s1','t0','t1','eNum'},fields))
    error('event struct does not contain the appropriate fields')
end


call = struct('s0',{},'s1',{},'t0',{},'t1',{},'cNum',{},'eNum',{});


%% iterate over each event
for eNum = 1:numel(event)
    E = event(eNum);
    
    % find starting sample for each block (if more than one)
    blk1 = (E.s0(1) : BLOCKSIZE : E.s1(1));
    blk2 = (E.s0(2) : BLOCKSIZE : E.s1(2));
    
    % split long events into multiple overlapping event blocks if block size exceeded
    [block(1:numel(blk1))] = deal(struct('s0',[0 0],'s1',[0 0]));
    for bNum = 1:numel(blk1)
        block(bNum).s0(1) = blk1(bNum);
        block(bNum).s1(1) = min(E.s1(1), blk1(bNum) + BLOCKSIZE);
        block(bNum).s0(2) = blk2(bNum);
        block(bNum).s1(2) = min(E.s1(2), blk2(bNum) + BLOCKSIZE);
        if any(block(bNum).s1 >= E.s1)
            block(bNum+1:end) = [];     % remove duplicate blocks at end
            break
        end
    end

    if numel(block) > 1;
        fprintf('Splitting event number %d into %d separate overlapping events\n',eNum,numel(block))
        fprintf('Processing block number    ')
    end


    %% iterate over each separated block
    for bNum = 1:numel(block)
        if numel(block) > 1
            fprintf('%s%.3d',char([8 8 8]),bNum);
        end
        
        B = block(bNum);
        idx1 = (B.s0(1):B.s1(1));
        idx2 = (B.s0(2):B.s1(2));
        
        % index data channels in sequential order for each board
        if isempty(array)
            chNum = 1:224;
            chIdx = chNum;
        else
            ch = array.ch + 112*(array.bd-1);   % get shuffle order
            [chNum,chIdx] = sort(ch);           % sort, but retain order
        end
        ch1 = chNum(chNum <= 112);
        ch2 = chNum(chNum > 112) - 112;
        
        % read block data samples from raw binary SRZ files
        [res1,hdr1] = read_SRZ(fname1,idx1,ch1);
        N1 = size(res1,1);
        assert(N1>0,'Side 1 returned no data!')

        res2 = read_SRZ(fname2,idx2,ch2);
        N2 = size(res2,1);
        assert(N2>0,'Side 2 returned no data!')

        % check to ensure both sides have equal length, truncate if necessary
        if (N1 ~= N2)
            if N1 > N2
                res1(N2+1:end,:) = [];
            else
                res2(N1+1:end,:) = [];
            end
            warning('AZ_SPLIT_EVENT:trunc','Data lengths unequal in block %d - truncating %d samples',bNum,abs(N1-N2))
        end
        
        res = [res1 res2];              % combine data from each board
        ts.data = zeros(size(res));     % init data struct
        ts.data(:,chIdx) = res;         % reverse channel sort order
        clear res*                      % free memory chunk
        
        % remove DC offset / high-pass filter data
        %ts.data = ts.data - ones(size(ts.data,1),1)*mean(ts.data);
        [b,a] = butter(2,10e3./hdr1.fs,'high');
        ts.data = filtfilt(b,a,ts.data);
        ts.fs = hdr1.fs;
        
        %% detect calls in data using instantaneous amplitude
        energy = abs(hilbert(ts.data));         % approximate with the magnitude of the hilbert transform
%        [x,y] = find(energy > gamma);           % search for threshold crossings
%        energy = sparse(x,y,1);                 % convert to sparse matrix
        energy = sum(energy,2);                 % sum energy over all channels
        eIdx = find(energy > gamma)';             % detect signals above threshold
        
        if 1
            figure(1)
            plot(energy)
            title(sprintf('Block %d of %d',bNum,numel(block)))
            drawnow
            
            figure(2)
            plot(ts.data)
            title('Timeseries data')
            drawnow
        end
        
        % if no detections found, continue to next block
        if isempty(eIdx)
            continue
        end
        
        % find start/stop samples (+/- padding) for one or multiple calls
        s0 = [eIdx(1) eIdx(1+find(diff(eIdx) > 472))] - nPad;
        s1 = [eIdx(find(diff(eIdx) > 472)) eIdx(end)] + nPad;
        
        % ensure start/stop is within block indices
        s0 = unique(max(s0,1));
        s1 = unique(min(s1,BLOCKSIZE));
        
        % append call struct with detected call(s)
        %[call(end+1:end+numel(s0))] = deal(struct('s0',[0 0],'s1',[0 0],'t0',[0 0],'t1',[0 0]));
        for c=1:numel(s0)
            % assign sample numbers
            C.s0 = B.s0 + s0(c) - 1;
            C.s1 = B.s0 + s1(c) - 1;
            
            % assign times using time from side 1
            C.t0 = E.t0 + (C.s0(1) ./ ts.fs);
            C.t1 = E.t0 + (C.s1(1) ./ ts.fs);
            
            % assign call/event number
            cNum = cNum + 1;
            C.cNum = cNum;
            C.eNum = eNum;
            
            % append call to data struct
            if numel(call) > 0 && any(call(end).s1 >= C.s0 - 1)     % detect if part of previous call
                cNum = cNum - 1;    % decrement counter
                call(end).s1 = C.s1;
                call(end).t1 = C.t1;
            else
                call(end+1) = C;
            end
        end
    end
    
    % remove last block if caused by noise
    if numel(call) && any(call(end).s1 >= E.s1)
        call(end) = [];
        cNum = cNum - 1;
    end
    
    %fprintf('   Detected %d calls in event %d\n', numel(c), eNum)

    % clear memory after each event is processed
    clear block
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save results to file
if ~isempty(callfile)
    fprintf('\nSaving event index to "%s"...\n',callfile)
    save(callfile,'call');
end
