function az_process_data(varargin)
% AZ_PROCESS_DATA  processes beam data and generates audio/video files for
% each trial found in the raw SRZ data sets
%
% az_process_data - prompts user for a folder.  All SRZ files located in
% the directory or subdirectories will be processed if a matching side1/2 
% pair is found

DEBUG = true;

% prompt user for working directory
wDir = uigetdir('Select a data directory');

% search directory for srz files
fnames = findfiles(wDir,'\.srz$');

% find which files are side1 by name
res = strfind(fnames,'side1.srz');
idx1 = find(~cellfun(@isempty,res));

%% iterate over each pair of files (side 1 and 2)
for n = 1:numel(idx1)
    % find prefix of filename
    [pname, fname, ext] = fileparts(fnames{idx1(n)});
    sLoc = regexpi(fname,'\_');
    prefix = fname(1:sLoc(end)-1);
    
    % verify side2 exists
    fname1 = fullfile(pname, [prefix '_side1' ext]);
    fname2 = fullfile(pname, [prefix '_side2' ext]);
    if ~existfile(fname2)
        warning('AZ_PROCESS_DATA:File Not Found','No matching side2 found for "%s" in "%s"\n  *** Skipping file ***',fullfile(fname,ext),pname)
        continue
    end
    
    %% load call map, if exists, otherwise create
    callfile = fullfile(pname, sprintf('%s_callmap.mat',prefix));
    if existfile(callfile)
        load(callfile,'callmap');
    else
        fprintf('Callmap not found.  Restarting ')
        callmap = az_detect(fname1,fname2);
        fprintf('\nSaving call index to "%s"...\n\n',callfile)
        save(callfile,'callmap');
    end
    
    %% separate data set into trials from call map
    %%%% eventually, this should use trigger count with a switch toggled for each trial
    t = [callmap(:).t0];                    % extract starting times of all events
    t = find(diff(t(1:2:end)) > 0.6) + 2;   % use difference between events as trial indicator
    t = t(diff(t) > 50);                    % remove trials with less than 50 pulses
    t = [1 t numel(callmap)+1];             % add first and last events

    %% iterate over each trial
    tic
    for n = 1:numel(t)-1
        try
            % show time series for each call
            if DEBUG
                plotTimeSeries(fname1,fname2,callmap(t(n):t(n+1)-1))
                drawnow
            end

            % process beam data
            beamfile = fullfile(pname, sprintf('%s_beams_%d-%d',prefix,t(n),t(n+1)-1));
            if ~existfile(beamfile)
                [beam, ref] = az_process_beams(fname1,fname2,t(n):t(n+1)-1,beamfile);
            else
                load(beamfile)
            end

            % generate audio track
            wavname = fullfile(pname, sprintf('%s_call_%d-%d.wav',prefix,t(n),t(n+1)-1));
            if ~existfile(wavname)
                genAudioTrack(ref, wavname);
            end

            % generate video track
            aviname = fullfile(pname, sprintf('%s_call_%d-%d.avi',prefix,t(n),t(n+1)-1)) ;
            if ~existfile(aviname)
                genVideoTrack(beam, ref, aviname);
            end

        catch
            % bail out and continue with next call
            warning('Holy Crap!')
            ERR = lasterror;
            disp(ERR.message)
            continue
        end
        
        % clear memory for next run
        clear beam ref
        
    end
    toc
    
end