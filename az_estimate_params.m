function ref = az_estimate_params(ref,ts,cmap)
% AZ_ESTIMATE_PARAMS  takes processed data from each call and derives basic
% bulk parameters, such as start time, pulse length, stop frequency, etc.

fprintf('\n***********************************************\n')
fprintf('Extracting call parameters from time series data\n')
ref.cNum = cmap.cNum;
ref.t0 = ts.t0 + cmap.t0(1);            % add absolute trigger time to relative start time
ref.t1 = ts.t1 + cmap.t1(1);            % add absolute trigger time to relative start time
ref.tlen = ts.tlen;
ref.ch = ts.refch;
ref.data = real(ts.data(:,ts.refch));     % save reference time series data
ref.fs = ts.fs;



%fprintf('Extracting call parameters from frequency domain data\n')



%fprintf('Extracting call parameters from processed beam data\n')
