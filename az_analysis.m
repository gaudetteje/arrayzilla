function fd = az_analysis(ts,varargin)
% AZ_ANALYSIS  Computes frequency response of each channel

% optional input parameters
MODE = 'ft';
if nargin > 1
    MODE = varargin{1};
end

fprintf('\n\n*****************************************\n')

% realign data samples given source position and time of arrival
%TBD

% estimate magnitude/phase per frequency bin
fprintf('Analyzing spectral content of received data\n')
switch lower(MODE)
    
    %%% mode 1 - use FT of each call
    case 'ft'
        
        % zero pad data
        ts.data = [zeros(250,size(ts.data,2)); ts.data; zeros(250,size(ts.data,2))];
        
        nfft = 2*floor(size(ts.data,1)/2);       % ensure even to avoid warning msg in calc_spectrum
        
        fd = calc_spectrum(ts,nfft,@blackmanharris,1,0.5,'half','ac','Vrms');

    %%% mode 2 - Hilbert analysis of each call
    case 'ht'
        % estimate the intstantaneous phase law for the call
        %IF = ;
        
        % remove echoes and reverb using pre-warping filter
        ts2 = mca_iffilt(ts.data,IF,ts.fs); %%% if this doesn't work on multiple signals then modify mca_iffilt to do so!
        
        % estimate parameters of remaining components
        [IMF, IA, IF] = mca_extract(ts2, ts.fs);
        
        % assign to frequency domain structure
        fd.freq = IF;
        fd.mag = IA(IF > .3*max(IF));
        fd.magdb = 20*log10(fd.mag);
    otherwise
        error('Unknown analysis mode used')
        
end
