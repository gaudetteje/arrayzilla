function ref = az_estimate_params(ts,event)
% AZ_ESTIMATE_PARAMS  takes processed data from each call and derives basic
% bulk parameters, such as start time, pulse length, stop frequency, etc.

fprintf('\n***********************************************\n')
fprintf('Extracting call parameters from time series data\n')
ref.eNum = event.eNum;
ref.t0 = ts.t0 + event.t0(1);            % add absolute trigger time to relative start time
ref.t1 = ts.t1 + event.t1(1);            % add absolute trigger time to relative start time
ref.tlen = ts.tlen;
ref.ch = ts.refch;
ref.data = real(ts.data(:,ts.refch));     % save reference time series data
ref.fs = ts.fs;


%fprintf('Extracting call parameters from frequency domain data\n')


%fprintf('Extracting call parameters from processed beam data\n')


ref.done = [];
ref.error = [];
