function ts = az_armaloss(ts, range)
% AZ_ARMALOSS  correct for transmission losses re 1m

PLOTFLAG = false;

Nch = size(ts.data,2);

% data params
fInt = 1e3;             % frequency domain sampling interval [Hz]

% system model params
Nb = 10;                % number of zeros
Na = 10;                % number of poles
maxIter = 10;           % maximum iterations with Steiglitz-McBride algorithm

% environmental params - should use measured parameters for more accurate results
h_r = 20;               % relative humidity [%]
T = 20;                 % ambient temperature [deg C]
P = 101.325;            % barometric pressure [kPa]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate frequency dependent transmission loss in air
f = (0:fInt:ts.fs/2);
alpha = calcAbsorptionCoef(f, h_r, T, P);

% calculate transmission losses for each channel
TL = zeros(length(f), Nch);
for cNum = 1:Nch
    TL(:,cNum) = 20*log10(range(cNum)) + alpha.*range(cNum);     % assume spherical (free-field) propagation
end

% define desired magnitude response from calculated transmission loss
TL = 10.^(TL/20);
TL = [TL; flipud(TL(2:end-1,:))];      % mirror full spectrum (removing point at Nyquist)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iterate over each channel
B = zeros(Nb+1,Nch);
A = zeros(Na+1,Nch);

fprintf('\n\n*****************************************\n')
fprintf('Inverting transmission losses on channel    ')
for cNum = 1:Nch
    fprintf('%s%.3d',char([8 8 8]),cNum);
    
    
    Hd = TL(:,cNum);            % get response of current channel
    
    % estimate system transfer function poles/zeros
    Hd = abs(Hd);               % ensure magnitude only data
    L = size(Hd,1);
    W = dftmtx(L); Wb = W; Wa = W;
    Wb(:,Nb+2:L) = []; Wa(:,Na+2:L) = [];

    % generate the autocorrelation function
    r = ifft(Hd.^2);

    % construct an initial system model by Levinson-Durbin (AR), follow with Prony (ARMA)
    aL = levinson(r,floor(L/2));
    hL = impz(1,aL,Nb+2*Na+2);
    [b,a] = prony(hL,Nb,Na);

    % iteratively refine pole/zero positions with frequency domain Steiglitz-McBride (ARMA)
    for i = 1:maxIter,
        [Hi,w] = freqz(b,a,L,'whole');
        Hai = freqz(1,a,L,'whole');
        Pi = exp(1i*angle(Hi));
        HdPi = Hd.*Pi;
        b = (diag(Hai)*Wb)\HdPi;
        a = (diag(HdPi.*Hai)*Wa)\(diag(Hai)*Wb*b);
    end

    % force real filter coefficients (Hd should be symmetric)
    if (sum(imag(b)) + sum(imag(a)) > 1e-10)
        warning('Poles and/or zeros not entirely real.  Possibly throwing away significant imaginary parts.')
        fprintf('Note:  The desired magnitude response, Hd, should be kept symmetric to ensure real coefficients!\n')
    end
    b = real(b.');
    a = real(a.');

    % scale all coefficients evenly, forcing a0=1
    b = b/a(1);
    a = a/a(1);
    
    % filter data through ARMA model
    ts.data(:,cNum) = filtfilt(b,a,ts.data(:,cNum));

end
fprintf('... Done!\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot frequency dependent absorption function
if PLOTFLAG
    calcAbsorptionCoef(f, h_r, T, P);
    
    figure
    [H,w4] = freqz(B,A,4*L,'whole');        % interpolated frequency response
    plot(w4/pi,abs(H),'k-',w/pi,TL,'k--')
    grid on;
    title('ARMA System Model, H(z) = S(z) = B(z)/A(z)')
    xlabel('Normalized frequency (w/pi)')
    ylabel('Magnitude')
    legend('Model magnitude','Desired magnitude','location','Best')

    H = H(1:4:4*L);     % downsample frequency domain

    figure
    plot(w/pi,Hd - abs(H),'r')
    grid on;
    title('ARMA System Model Error, Hd(z) - H(d)')
    xlabel('Normalized frequency (w/pi)')
    ylabel('Magnitude')

    figure;
    zplane(b,a);
    title('Pole/zero plot')

    fprintf('Sum of sq. error:  %g\n', sum((Hd - abs(H)).^2))
end