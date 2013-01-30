function [ts, call] = az_split_event(fname1,fname2,event,array)
% AZ_SPLIT_EVENT searches raw binary data for one or more calls in each
%     event
%
% CALLS = az_convert(FNAME1,FNAME2,EVENT) retrieves the data in the
%     SRZ file pair (FNAME1, FNAME2) using the CALLMAP structure.
%
% The function returns a new Mx1 structure, CALLS, with identical fields
%     as EVENT.
%
% Note:  EVENT is generated using az_detect and ARRAY is generated using
% az_positions and az_channelmap.

fprintf('\n***********************************************\n')

% set default parameters
BLOCKSIZE = 23666;              % process blocks of 100ms maximum
OVERLAP = floor(BLOCKSIZE*0.1); % use 10% overlap to avoid cutting off calls
gamma = 1;

%% split long events into overlapping blocks

% find starting sample for each block (if more than one)
blk1 = (event.s0(1) : BLOCKSIZE-OVERLAP : event.s1(1));
blk2 = (event.s0(2) : BLOCKSIZE-OVERLAP : event.s1(2));

% split 1st event into multiple overlapping events if block size exceeded
[call(1:numel(blk1))] = deal(struct('s0',[0 0],'s1',[0 0],'t0',[0 0],'t1',[0 0]));
for e = 1:numel(blk1)
    
    call.s1(1) = min(event.s1(1), blk1 + BLOCKSIZE);
    call.s1(2) = min(event.s1(2), blk2 + BLOCKSIZE);
end


%% read first event data from raw binary SRZ files
% find sample numbers to retrieve on each board
idx1 = (call.s0(1):call.s1(1));
idx2 = (call.s0(2):call.s1(2));

% index data channels in sequential order for each board
chNum = array.ch + 112*(array.bd-1); % get shuffle order
%%% try doing this using array.ch
%%% sort, read_SRZ, unsort
ch1 = 1:112;
ch2 = 1:112;

% read in both data sets
res1 = read_SRZ(fname1,idx1,ch1);
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
    warning('AZ_CONVERT:trunc','Data lengths unequal - truncating %d samples',abs(N1-N2))
end

res = [res1 res2];                  % combine data from each board
ts.data = zeros(size(res,1), length(array.ch));         % init data struct

%%%%
% reshuffle 
for n=1:length(chNum)
    %fprintf('ts.data(:,%d) = res(:,%d) [Ch %d, Bd %d]\n',n,chIdx(n),array.ch(n),array.bd(n))
    ts.data(:,n) = res(:,chNum(n));
end
%%%%

% remove DC offset
ts.data = ts.data - ones(size(ts.data,1),1)*mean(ts.data);
ts.fs = hdr1.fs;


%% detect calls in data
energy = ts.data.^2;                         % convert to energy

% detect calls using energy
valid = ts.data > gamma;
% prepend events struct with valid call(s)
%TBD

% reduce time series to valid data sections
%TBD

% throw warning if no valid calls were detected in block
if numel(call) == 0
    warning('AZ_CONVERT:detect','No valid calls were detected')
end

%% restructure event list to include each detected call
%%% should return index of each block event to avoid reading everything in
% read first block, return ts and remaining blocks (if any)
