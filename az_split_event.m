function [call, ts] = az_split_event(fname1,fname2,event,varargin)
% AZ_SPLIT_EVENT searches raw binary data for one or more calls in each event
%
% CALL = az_split_event(FNAME1,FNAME2,EVENT) returns a CALL index struct
%       from data in the SRZ file pair (FNAME1, FNAME2) referenced in
%       EVENT struct
% CALL = az_split_event(FNAME1,FNAME2,EVENT,ARRAY) uses the channel
%       mapping in ARRAY to eliminate dead channels from analysis
% CALL = az_split_event(...,CALLFILE) entering a string as the last
%       parameter will save the resulting call struct to CALLFILE
%
% The function returns a new Mx1 structure, CALLS, with identical fields
%     as EVENT.  Long events (100ms or greater) are split into multiple
%     overlapping events before searching for calls.  This is to avoid
%     memory issues and allows processing continuous recorded data where no
%     trigger events are found.
%
% Note:  EVENT is generated using az_detect and ARRAY is generated using
% az_positions and az_channelmap.

% TODO:
% This is the first time we actually read data from the array channels
% - detect bad channels using sum(E,1) to look across time
% - iterate over each data block to detect ALL calls, not just first block
% - avoid reading all 224 channels - reorder and read subset of good
%   channels (TIME THIS TO ENSURE EFFICIENCY)

fprintf('\n***********************************************\n')

% set default parameters
BLOCKSIZE = 23666;              % process blocks of 100ms maximum
OVERLAP = floor(BLOCKSIZE*0.1); % use overlap to avoid cutting off calls

gamma = 1e-2;                   % normalized threshold for amplitude detection
nCh = 5;                        % number of channels required for threshold
nPad = 200;                     % number of samples to pad around detected calls
DEBUG = false;

% init parameters
array = [];
cNum = 0;                       % detected call counter
callfile = [];

% handle optional inputs
switch nargin
    case 5
        array = varargin{1};
        callfile = varargin{2};
    case 4
        res = varargin{1};
        if ischar(res)
            callfile = res;
        else
            array = res;
        end
    case 3
    otherwise
        error('Incorrect number of parameters entered')
end

call = struct('s0',{},'s1',{},'t0',{},'t1',{},'cNum',{},'eNum',{}); % init empty struct (populated later)
assert(OVERLAP/BLOCKSIZE <= 0.5,'OVERLAP should be between 0% and 50%')


%% split long events into overlapping blocks
for eNum = 1:numel(event)
    E = event(eNum);
    
    % find starting sample for each block (if more than one)
    blk1 = (E.s0(1) : BLOCKSIZE-OVERLAP : E.s1(1));
    blk2 = (E.s0(2) : BLOCKSIZE-OVERLAP : E.s1(2));

    % split each event into multiple overlapping event blocks if block size exceeded
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
    end


    %% iterate over each separated block
    for bNum = 1:numel(block)
        
        B = block(bNum);
        idx1 = (B.s0(1):B.s1(1));
        idx2 = (B.s0(2):B.s1(2));

        % index data channels in sequential order for each board
        if isempty(array)
            chNum = 1:224;
        else
            %%% try doing this using array.ch and array.bd values
            %%% i.e. sort, read_SRZ, unsort
            chIdx = array.ch + 112*(array.bd-1);    % get shuffle order
            %[chNum,chIdx2] = sort(chIdx); %% NEED TO VERIFY SHUFFLE ORDER
            %GETS MAINTAINED SO WE CAN REVERSE IT!
        end
        ch1 = 1:112;
        ch2 = 1:112;

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

        res = [res1 res2];                  % combine data from each board
        ts.data = zeros(size(res,1), length(chNum));         % init data struct

        %%%%
        % reshuffle 
        for n=1:length(chNum)
            %fprintf('ts.data(:,%d) = res(:,%d) [Ch %d, Bd %d]\n',n,chIdx(n),array.ch(n),array.bd(n))
            ts.data(:,n) = res(:,chNum(n));
        end

        % remove DC offset
        ts.data = ts.data - ones(size(ts.data,1),1)*mean(ts.data);
        ts.fs = hdr1.fs;

        %% detect calls in data using energy
        %E = ts.data.^2 * 2^-32;         % convert to normalized energy
        energy = abs(ts.data) * 2^-15;       % normalize abs val
        [x,y] = find(energy > gamma);        % search for threshold crossings
        energy = sparse(x,y,1);              % convert to sparse matrix
        energy = sum(energy,2);                   % sum threshold crossings over all channels
        eIdx = find(energy > nCh);               % detect calls if 10 or more signals coincide

        if DEBUG
            fprintf('.')
            figure(1)
            plot(energy)
            title(sprintf('Block %d of %d',bNum,numel(block)))
            drawnow
        end
        
        % if no detections found, continue to next block
        if isempty(eIdx)
            continue
        end

        % find start/stop samples (+/- padding) for one or multiple calls
        s0 = [eIdx(1) eIdx(1+find(diff(eIdx) > 472))] - nPad;
        s1 = [eIdx(find(diff(eIdx) > 472)) eIdx(end)] + nPad;

        % if call overlaps the start or end of a block, ignore and move to next overlapping block
        if any((s0 < 1) & (bNum > 1))
            warning('AZ_SPLIT_EVENT:startblock','call starts before block %d',bNum)
            continue
        end
        if any(s1 > size(energy,1))
            warning('AZ_SPLIT_EVENT:endblock','call ends after block %d',bNum)
            continue
        end
        %%% TBD - REMOVE REDUNDANT CALLS HERE %%%

        % append call struct with detected call(s)
        %[call(end+1:end+numel(s0))] = deal(struct('s0',[0 0],'s1',[0 0],'t0',[0 0],'t1',[0 0]));
        for c=1:numel(s0)
            cNum = cNum + 1;
            C.s0 = B.s0 + s0(c) - 1;
            C.s1 = B.s0 + s1(c) - 1;
            C.t0 = E.t0 + (C.s0 ./ ts.fs);
            C.t1 = E.t0 + (C.s1 ./ ts.fs);
            C.cNum = cNum;
            C.eNum = eNum;
            
            % append call to data struct
            call(end+1) = C;
        end
    end
    fprintf('   Found %d calls in event %d\n',numel(call),eNum)

    % clear memory after each event is processed
    clear block
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save results to file
if ~isempty(callfile)
    fprintf('\nSaving event index to "%s"...\n',callfile)
    save(callfile,'call');
end
