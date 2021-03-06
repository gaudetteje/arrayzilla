function fd = az_analysis(ts,varargin)
% AZ_ANALYSIS  Computes frequency response of each channel

fprintf('\n***********************************************\n')

% estimate magnitude/phase per frequency bin
fprintf('Analyzing spectral content of received data\n')

% FFT analysis on zero-padded data
x.fs = ts.fs;
x.data = [zeros(250,size(ts.data,2)); ts.data; zeros(250,size(ts.data,2))];
nfft = 2*floor(size(x.data,1)/2);       % ensure even to avoid warning msg in calc_spectrum
fd = calc_spectrum(x,nfft,@blackmanharris,1,0.5,'half','ac','Vrms');
%%% Try using peak value of STFT since there is an inherent smoothing of
%%% the TF plane...  This should average out some of the noise in each
%%% channel

% Hilbert analysis of each call
if isfield(ts,'fm1')
    % initialize matrices
    IA1 = zeros(size(ts.fm1));
    IA2 = zeros(size(ts.fm1));
    IF1 = zeros(size(ts.fm1));
    IF2 = zeros(size(ts.fm1));

    % estimate the instantaneous phase law for harmonics in each call
    for n=1:size(ts.fm1,2)
        [~,IA1(:,n),IF1(:,n)] = mca_extract(ts.fm1(:,n), ts.fs, 0);
        [~,IA2(:,n),IF2(:,n)] = mca_extract(ts.fm2(:,n), ts.fs, 0);
    end

    % assign to frequency domain structure
    fd.if1 = IF1;
    fd.ia1 = IA1;

    fd.if2 = IF2;
    fd.ia2 = IA2;
end
