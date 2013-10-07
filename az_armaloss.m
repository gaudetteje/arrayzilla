function ts = az_armaloss(ts, range, varargin)
% AZ_ARMALOSS  correct for transmission losses re 1m

PLOTFLAG = false;
if nargin > 2
    PLOTFLAG = varargin{1};
end

Nch = size(ts.data,2);

% data params
fInt = 1e3;             % frequency domain sampling interval [Hz]

r0 = 1;                 % SPL reference distance [m]

% environmental params - should use measured parameters for more accurate results
h_r = 20;               % relative humidity [%]
T = 20;                 % ambient temperature [deg C]
P = 101.325;            % barometric pressure [kPa]

% system model params
Nb = 8;                % number of zeros
Na = 8;                % number of poles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate frequency dependent transmission loss in air
f = (0:fInt:ts.fs/2);
alpha = calcAbsorptionCoef(f, h_r, T, P);

% calculate transmission losses for each channel
TL = zeros(length(f), Nch);
for cNum = 1:Nch
    TL(:,cNum) = 20*log10(range(cNum)./r0) + alpha.*(range(cNum)-r0);     % assume spherical (free-field) propagation
end

% define desired magnitude response from calculated transmission loss
TL = 10.^(TL/20);                       % define inverse magnitude response in linear units
TL = [TL; flipud(TL(2:end-1,:))];       % mirror spectrum
TL = sqrt(TL);                          % take square root since forward-backward filtering squares magnitude response

%% plot frequency dependent absorption function
if PLOTFLAG
    calcAbsorptionCoef(f, h_r, T, P);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iterate over each channel
fprintf('\n***********************************************\n')
fprintf('Inverting transmission losses on channel    ')
for cNum = 1:Nch
    fprintf('%s%.3d',char([8 8 8]),cNum);
    
    % get double-sided frequency response of current channel
    Hd = TL(:,cNum);
    
    % use frequency-domain StMcB IIR design method
    [b,a] = armacmpx2(Hd,Nb,Na,true);
    
    % filter data through ARMA model
    ts.data(:,cNum) = filtfilt(b,a,ts.data(:,cNum));

end
fprintf('... Done!\n')
