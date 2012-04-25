function [call_ts]=az_findcalls(DATA)
% search for calls on a specified channel over a specified range of time
% 
%Settings
refperiod=.01;
refperiod=refperiod*Fs;


%% find calls
ma_windowsize=200;
data=filter_bandpass(Wf,15e3,70e3,Fs);
data=filtfilt(ones(1,ma_windowsize)/ma_windowsize,1,abs(data));
data=data/(median(data));


type='findpeaks';

switch lower(type)
    case {'loop'}
        refcount=refperiod;
        call_Ts=[];
        CALLSTATE=0;
        for k3=1:length(data)
            refcount=refcount+1;
            if CALLSTATE==0
                if (data(k3)>3 && refcount>refperiod)
                    ctime=k3/Fs;
                    call_Ts=[call_Ts ctime];
                    refcount=0;
                    CALLSTATE=1;
                end
            elseif CALLSTATE==1
                if data(k3)<1.5
                    CALLSTATE=0;
                end
            end
            
            
        end
        
    case {'findpeaks'}
        [pks call_Ts]=findpeaks(data,'minpeakheight',3,'minpeakdistance',refperiod);
        call_Ts=call_Ts/Fs;
end


call_Ts=call_Ts+STARTTIME;
if 1
    T=(1:length(Wf))/Fs;
    plot(T,data)
    hold on
    plot(call_Ts,zeros(size(call_Ts)),'.r')
end