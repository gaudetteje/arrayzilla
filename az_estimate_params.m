function ref = az_estimate_params(ts,fd,beam,event)
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

% iterate over each frequency in beam
N = numel(beam.f);
ref.mra.f = beam.f;             % beam frequency index
ref.mra.pk = nan(N,1);          % main response axis - peak value
ref.mra.az = nan(N,1);          % main response axis - azimuth
ref.mra.el = nan(N,1);          % main response axis - elevation
for f = 1:N
    [x,y] = find(beam.FFT(:,:,f) == max(max(beam.FFT(:,:,f))));
    ref.mra.pk(f) = beam.FFT(x,y,f);
    ref.mra.az(f) = beam.az(y);
    ref.mra.el(f) = beam.el(x);
end

% flags for debugging use
ref.done = [];
ref.error = [];
