function fd = az_analysis(ts,varargin)
% AZ_ANALYSIS  Computes frequency response of each channel

% optional input parameters
%MODE = 'ft';
%if nargin > 1
%    MODE = varargin{1};
%end

fprintf('\n\n*****************************************\n')

% estimate magnitude/phase per frequency bin
fprintf('Analyzing spectral content of received data\n')
%switch lower(MODE)
    
    %%% mode 1 - use FT of each call
    %case 'ft'
        
        % zero pad data
        ts.data = [zeros(250,size(ts.data,2)); ts.data; zeros(250,size(ts.data,2))];
        
        nfft = 2*floor(size(ts.data,1)/2);       % ensure even to avoid warning msg in calc_spectrum
        
        fd = calc_spectrum(ts,nfft,@blackmanharris,1,0.5,'half','ac','Vrms');
    
    %%% mode 2 - Hilbert analysis of each call
    %case 'ht'
    
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
        %fd.ia1db = db(fd(1).mag);
        
        fd.if2 = IF2;
        fd.ia2 = IA2;
        %fd.magdb = db(fd(2).mag);
    end
    
%    otherwise
%        error('Unknown analysis mode used')
    
%end
