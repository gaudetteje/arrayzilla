function ts = az_convert(fname1,fname2,call,array)
% AZ_CONVERT reads data from the raw binary files into a MATLAB structure
%
% TS = az_convert(FNAME1,FNAME2,CALLMAP,ARRAY) retrieves the data in the
% SRZ file pair (FNAME1, FNAME2) using the CALLMAP structure.  Channel
% numbers are mapped to array indices and coordinates using the ARRAY
% struct.
%
% The function returns a structure, TS, with fields 'data' and 'fs'.
% 'data' is an MxN matrix with M samples across N channels.  'fs' is a
% scalar value holding the data sampling rate.
%
% Note:  CALLMAP is generated using az_detect and ARRAY is generated using
% az_positions and az_channelmap.

fprintf('\n\n***********************************************\n')

% find sample numbers to retrieve on each board
idx1 = (call.s0(1):call.s1(1));
idx2 = (call.s0(2):call.s1(2));

% index data channels in sequential order for each board
ch1 = 1:112;
ch2 = 1:112;

% read in both data sets
fprintf('Reading call data from side 1...  \n')
[res1, hdr1] = read_SRZ(fname1,idx1,ch1);
N1 = size(res1,1);
assert(N1>0,'Side 1 returned no data!')
fprintf('  Done!\n\n')

fprintf('Reading call data from side 2...  \n')
res2 = read_SRZ(fname2,idx2,ch2);
N2 = size(res2,1);
assert(N2>0,'Side 2 returned no data!')
fprintf('  Done!\n\n')

% check to ensure both sides have equal length, truncate if necessary
fprintf('Rearranging data columns to match array definition\n')
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
chNum = array.ch + 112*(array.bd-1); % get shuffle order
for n=1:length(chNum)
    %fprintf('ts.data(:,%d) = res(:,%d) [Ch %d, Bd %d]\n',n,chIdx(n),array.ch(n),array.bd(n))
    ts.data(:,n) = res(:,chNum(n));
end

ts.data = ts.data * (5*2^-16);      % convert to Volts
ts.data = ts.data - 2.5;            % remove DC offset

ts.fs = hdr1.fs;
