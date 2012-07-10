function [coords stats]=TDOA_frame(WFMAT,MICPOS,Fs,varargin)
%This function uses the TDOA reconstruction algorithm proposed in Silvermen
%et al 2008 in order to reconstruct the position of the bat using n
%microphones.
%
% COORDS = TDOA_FRAME(WFMAT,MICPOS,FS)  takes a numCh x numSamp data matrix,
%          WFMAT, a numCh x 3 sensor position matrix, MICPOS, and sampling rate, FS
%

if ~isempty(varargin)
    REFCHANNEL=varargin{1};
    assert(REFCHANNEL <= size(MICPOS,1),'Reference channel exceeds number of receive channels')
end

%% Settings:
c=calcSoundSpeed(26);

%thresholds
ampthresh=.0005; %amplitude threshold to accept channel
corrthresh=.0005; %cross-correleation max threshold to accept channel
residthresh=.05; %resitudal threshold to accept position estimate

%othersettings
n_compchannels=6;

% filter call
filtertest=0;

%sample test
sampletest=0;
if sampletest==1
    FIG(1)=figure;
    FIG(2)=figure;
end

%Ploting options
plottest=0; % set to 1 to plot, or 0 to not plot
if plottest==1
    fig(1)=figure('position',[10 10 500 900]);
    fig(2)=figure('position',[650 10 800 900]);
end

%Determine microphone sets (for the selection of a reference mics)


%% Select and Load Data

nchannels=size(WFMAT,1);
MSET(1).set=[1:round(nchannels/2)];
%MSET(2).set=[round(nchannels/2):nchannels];

%% Preallocate
stats.micuserec=zeros(1,nchannels);


%% Grab Call

for kch=1:nchannels %loop through channels populating Wf matrix
    if isequal(filtertest,1)
        WFMAT(kch,:)=wf_call_filter(WFMAT(kch,:),Fs);
    end
    callamp(kch)=sqrt(mean(WFMAT(kch,:).^2)); %calculate RMS amplitude on each channel
end



%%  calculate cross correlations

for kmset=1:length(MSET)
    if exist('REFCHANNEL')
        MSET(kmset).rmic=REFCHANNEL(kmset);
    else
        [mx idx]=max(callamp(MSET(kmset).set));  %pick mic with greatest amplitude as reference from each mic set
        MSET(kmset).rmic=MSET(kmset).set(idx);
    end
    
    MSET(kmset).t_ref_detect=0;
    
    MSET(kmset).s_ref_detect=(MSET(kmset).t_ref_detect+.01)*Fs;
    
    % perform cross correlations to find time differences
    if isfield(MSET,'corr')
        MSET=rmfield(MSET,'corr');
    end
    
    for k=1:nchannels
        [MSET(kmset).tno(k,1) MSET(kmset).corr(k,:) MSET(kmset).lags MSET(kmset).corrmax(k,1)]=wfcrosscorr_fa(WFMAT(MSET(kmset).rmic,:),WFMAT(k,:),Fs);
    end
    
    MSET(kmset).s_rel_corr(:,1)=MSET(kmset).tno(:,1)*Fs+MSET(kmset).s_ref_detect; 
end

if plottest==1
    figure(fig(1))
    
    T=(1:length(WFMAT(1,:)))/Fs;
    for kch=1:nchannels
        
        
        plot(T*1e3,kch+8*WFMAT(kch,:),'k')
        %ylim([-.1 .1])
        hold on
        
        ylabel('Channel')
        xlabel('Time (ms)')
    end
    
    figure(fig(2))
    nset=length(MSET);
    for kmset=1:nset
        
        
        subplot(1,nset,kmset)
        plot(MSET(kmset).tno*1e3,1:nchannels,'.k')
        hold on
        plot(MSET(kmset).tno(MSET(kmset).rmic)*1e3,MSET(kmset).rmic,'.r')
        hold off
        
        ylabel('Channel')
        xlabel('Measured TDOA (ms)')
        xlim([-2 2])
    end
    
end


%% Perform reconstruction

micdata.locations=MICPOS;
micdata.locations=micdata.locations(1:nchannels,:);

n_rows=(nchannels-1)*length(MSET);
voidrows=[];

for kmset=1:length(MSET)
    voidrows=[voidrows ((kmset-1)*nchannels+MSET(kmset).rmic)];  %void rows are ignored in the computation. The reference mike must be ignored
    for k=1:nchannels
        if (callamp(k)<ampthresh || MSET(kmset).corrmax(k)<corrthresh)    %if amplitude or correlation max is lower than threshold, mark row as void
            voidrows=[voidrows ((kmset-1)*nchannels+k)];
        end
    end
end

n_rows=nchannels*length(MSET)-length(voidrows);
disp(['Starting with  ' num2str(n_rows) ' correlations from ' num2str(nchannels) ' channels'])
% Generate matrix

for kmset=1:length(MSET)
    MSET(kmset).dno=MSET(kmset).tno*c; %dno is a term
    
    MSET(kmset).xomxn=micdata.locations(MSET(kmset).rmic,1)-micdata.locations(1:nchannels,1); %calculate three terms
    MSET(kmset).yomyn=micdata.locations(MSET(kmset).rmic,2)-micdata.locations(1:nchannels,2);
    MSET(kmset).zomzn=micdata.locations(MSET(kmset).rmic,3)-micdata.locations(1:nchannels,3);
    
    MSET(kmset).wno=(MSET(kmset).dno.^2 - micdata.locations(:,1).^2 + micdata.locations(MSET(kmset).rmic,1).^2 ... %calculate wno
        -micdata.locations(:,2).^2 + micdata.locations(MSET(kmset).rmic,2).^2 ...
        -micdata.locations(:,3).^2 + micdata.locations(MSET(kmset).rmic,3).^2)/2;
end


Afull=zeros(nchannels*length(MSET),3+length(MSET)); %Preallocate Matrix
wnofull=zeros(nchannels*length(MSET),1);
for kmset=1:length(MSET)
    Afull((nchannels*(kmset-1)+1):nchannels*kmset,1)=MSET(kmset).xomxn; %load matrix
    Afull((nchannels*(kmset-1)+1):nchannels*kmset,2)=MSET(kmset).yomyn;
    Afull((nchannels*(kmset-1)+1):nchannels*kmset,3)=MSET(kmset).zomzn;
    Afull((nchannels*(kmset-1)+1):nchannels*kmset,3+kmset)=MSET(kmset).dno;
    
    
    wnofull((nchannels*(kmset-1)+1):nchannels*kmset,1)=(MSET(kmset).dno.^2 - micdata.locations(:,1).^2 + micdata.locations(MSET(kmset).rmic,1).^2 ... %calculate wno
        -micdata.locations(:,2).^2 + micdata.locations(MSET(kmset).rmic,2).^2 ...
        -micdata.locations(:,3).^2 + micdata.locations(MSET(kmset).rmic,3).^2)/2;
end



%% Sample
if sampletest==1
    totalrows=nchannels*length(MSET);
    for kgen=1:2000
        %select channels to reconstruct randomly
        goodrows=randchoose(setdiff(1:totalrows,voidrows),6);
        voidrows_temp=setdiff(1:nchannels,goodrows); %remove all but the good rows
        
        A=Afull;
        wno=wnofull;
        voidrows_temp=sort(voidrows_temp,'descend');  %remove void rows
        for kdel=1:length(voidrows_temp)
            A(voidrows_temp(kdel),:)=[];
            wno(voidrows_temp(kdel),:)=[];
        end
        
        Apsi=pinv(A);  %calculate pseudoinverse
        X=Apsi*wno;   %calculate location
        X=X';
        X_sample(kgen,:)=X;
    end
    
    figure(FIG(1));
    subplot(3,1,1)
    hist(X_sample(:,1),100)
    subplot(3,1,2)
    hist(X_sample(:,2),100)
    subplot(3,1,3)
    hist(X_sample(:,3),100)
    
    figure(FIG(2));
    hold off
    
    plot3(X_sample(:,1),X_sample(:,2),X_sample(:,3),'.k')
    
    
    
end




%% Climb
if n_rows > n_compchannels %if more than 6 rows remain non-void, perform reconstruction
    
    %%%% Using Residual Find the best 6 channels from which to calculate the position %%%
    
    while n_rows > n_compchannels
        
        good_rows=setdiff(1:(nchannels)*length(MSET),voidrows);
        residual_test=zeros(1,length(good_rows));
        
        for krow=1:length(good_rows) %loop through each of the good rows
            
            A=Afull;
            wno=wnofull;
            
            
            voidrows_temp=[voidrows good_rows(krow)]; % temporary list of voidrows including the test row (krow)
            
            voidrows_temp=sort(voidrows_temp,'descend');  %remove void rows
            for kdel=1:length(voidrows_temp)
                A(voidrows_temp(kdel),:)=[];
                wno(voidrows_temp(kdel),:)=[];
            end
            
            Apsi=pinv(A);  %calculate pseudoinverse
            X=Apsi*wno;   %calculate location
            X=X';
            
            %calculate residual
            residualmatrix_temp=zeros(nchannels*length(MSET),1);
            dno_recalc=zeros(nchannels*length(MSET),1);
            
            for kmset=1:length(MSET)
                
                for k=1:nchannels
                    NUM=(kmset-1)*nchannels+k;
                    if ~ismember(NUM,voidrows_temp)
                        dno_recalc(NUM,1)=norm(X(1,1:3)-micdata.locations(k,1:3))-norm(X(1,1:3)-micdata.locations(MSET(kmset).rmic,1:3));
                        residualmatrix_temp(NUM,1)=abs(dno_recalc(NUM,1)-MSET(kmset).dno(k,1));
                    end
                end
            end
            
            % Save the residual without krow
            residual_test(krow)=sum(residualmatrix_temp);
            
        end
        
        % find and remove the worst row
        %disp(length(residual_test))
        %disp(length(good_rows))
        if 0
            
            figure(fig(3))
            subplot(2,1,1)
            plot(residual_test)
            title(['call ' num2str(kabscall) ' ncompchannels ' num2str(n_rows)])
            subplot(2,1,2)
            
            disp(X)
        end
        
        [mn idx]=min(residual_test);
        worstrow=good_rows(idx);
        
        voidrows=[voidrows worstrow];
        
        n_rows=nchannels*length(MSET)-length(voidrows);
        
    end
    
    %%% calclulate the actual position %%%
    good_rows=setdiff(1:nchannels,voidrows);
    
    
    A=Afull;
    wno=wnofull;
    voidrows=sort(voidrows,'descend');  %remove void rows
    for kdel=1:length(voidrows)
        A(voidrows(kdel),:)=[];
        wno(voidrows(kdel),:)=[];
    end
    
    
    Apsi=pinv(A);  %calculate pseudoinverse
    X=Apsi*wno;    %calculate location
    X=X';
    
    if 0
        figure(fig(3))
        plot(A(:,4),'.k')
        hold on
        plot(A(:,5),'.k')
        hold off
    end
    
    %save the residual
    %calculate residual
    residuaL.matrix=zeros(nchannels*length(MSET),1);
    dno_recalc=zeros(nchannels*length(MSET),1);
    
    for kmset=1:length(MSET)
        
        for k=1:nchannels
            NUM=(kmset-1)*nchannels+k;
            if ~ismember(NUM,voidrows_temp)
                dno_recalc(NUM,1)=norm(X(1,1:3)-micdata.locations(k,1:3))-norm(X(1,1:3)-micdata.locations(MSET(kmset).rmic,1:3));
                residual.matrix(NUM,1)=abs(dno_recalc(NUM,1)-MSET(kmset).dno(k,1));
            end
        end
    end
    
    
    
    
    idxs=find(residual.matrix);
    residual.sum=sum(residual.matrix(idxs));
    residual.mics=idxs;
    residual.mean=mean(residual.matrix(idxs));
    
    stats.micuserec=stats.micuserec+ismember(1:nchannels,residual.mics);
    
%    disp(['resid sum @ ' num2str(n_compchannels) ' = ' num2str(residual.sum)])
%    disp(['resid mean = ' num2str(residual.mean)])
%    disp(X)
%    disp([' '])

    
    
    % decide whether to keep call
    coords=X; %save coordinates
    if (residual.mean >= residthresh)
        warning('Mean of residuals is greater than threshold - convergence may have failed')
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SubFunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dtime,c,lags,corrmax]=wfcrosscorr_fa(wf1,wf2,Fs)


% Correlation signal test


[c,lags] = xcov(wf2,wf1);
[corrmax index]=max(c);
%[PKS LOC]=findpeaks(c,'minpeakheight',.1*corrmax,'minpeakdistance',.1*Fs);
%plot(lags,c)
%hold on
%plot(lags(LOC),PKS,'.r')


dsamplesrec=lags(index);

dtime=dsamplesrec/Fs;



