function ts = az_convert(fname1,fname2,call,array)
% AZ_CONVERT takes the raw binary files and converts to Voltage units

fprintf('\n\n***********************************************\n')

% find sample numbers to retrieve on each board
idx1 = (call.s0:call.s1);
idx2 = idx1;

% index data channels in sequential order for each board
ch1 = 1:112;
ch2 = 1:112;

% read in both data sets
fprintf('Reading call data from side 1...  \n')
[res1, hdr1] = read_SRZ(fname1,idx1,ch1);
assert(numel(res1)>0,'Side 1 returned no data!')
fprintf('  Done!\n\n')

fprintf('Reading call data from side 1...  \n')
[res2] = read_SRZ(fname2,idx2,ch2);
assert(numel(res2)>0,'Side 2 returned no data!')
fprintf('  Done!\n\n')

fprintf('Rearranging data columns to match array definition\n')
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
