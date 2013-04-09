function az_process_data(varargin)
% AZ_PROCESS_DATA  processes beam data and generates audio/video files
%
% az_process_data - prompts user for a folder.  All SRZ files located in
%   the directory or subdirectories will be processed if a matching side1
%   and side2 pair is found
%
% az_process_data(PATHNAME) searches through the path specified


% start logging text session
logname = [datestr(now,'yyyymmdd_HHMM') '_log.txt'];
diary(logname)
disp(repmat('#',1,70))
fprintf('Started processing data files:  %s\n',datestr(now))


DEBUG = false;

plotFreq = [30e3 60e3 90e3];       % beam frequencies to plot


%% prompt user for working directory
if ~nargin
    wDir = uigetdir('Select a data directory');
else
    wDir = varargin{1};
end

% search directory for srz files
fnames = findfiles(wDir,'\.srz$');

% find which files are side1 by name
res = strfind(fnames,'side1.srz');
idx1 = find(~cellfun(@isempty,res));

if isempty(idx1)
    fprintf('No data files found.  Aborting...\n\n')
end


%% iterate over each pair of files (side 1 and 2)
for m = 1:numel(idx1)
    try
        tic
        
        % find prefix of filename
        [pname, fname, ext] = fileparts(fnames{idx1(m)});
        sLoc = regexpi(fname,'\_');
        prefix = fname(1:sLoc(end)-1);

        % verify side2 exists
        fname1 = fullfile(pname, [prefix '_side1' ext]);
        fname2 = fullfile(pname, [prefix '_side2' ext]);
        
        disp(repmat('#',1,70))
        fprintf('Processing data files with prefix "%s":\n\t%s\t%s\n\n',prefix,fname1,fname2)

        if ~existfile(fname2)
            warning('AZ_PROCESS_DATA:File Not Found','No matching side2 found for "%s" in "%s"\n  *** Skipping file ***',fullfile(fname,ext),pname)
            continue
        end
        
        
        %% validate and realign data files as necessary
        disp(repmat('#',1,70))
        if isempty(findfiles('.','event\.mat$'))
            fprintf('Verifying and realigning SRZ data files\n')
            az_align_data(fname1,'auto');
            az_align_data(fname2,'auto');
        end
        
        %% load array structure, if exists, otherwise create it
        disp(repmat('#',1,70))
        arrayfile = fullfile(pname, [prefix '_array.mat']);
        if existfile(arrayfile)
            fprintf('Loading precomputed array definition...\n')
            load(arrayfile,'array')
        else
            fprintf('Defining array structure\n')
            array = az_define_array(arrayfile);
        end
        
        
        %% load event map, if exists, otherwise create it
        disp(repmat('#',1,70))
        eventfile = fullfile(pname, [prefix '_event.mat']);
        hdrfile = fullfile(pname, [prefix '_hdr.mat']);
        if ~existfile(eventfile)
            fprintf('Detecting events...\n')
            event = az_detect_events(fname1,fname2,eventfile,hdrfile);
        end

        callfile = fullfile(pname, [prefix '_call.mat']);
        if existfile(callfile)
            fprintf('Loading precomputed call structure...\n')
            load(callfile,'call');
        else
            fprintf('Loading precomputed event structure...\n')
            load(eventfile,'event');
            fprintf('Call map not found.  Parsing data now...')
            call = az_split_event(fname1,fname2,event,array,callfile);
        end
        
        % show time series for each event
        if DEBUG
            fprintf('Plotting time series of all calls\n')
            plotTimeSeries(fname1,fname2,call)
            drawnow
        end

        % process beam data
        disp(repmat('#',1,70))
        beamfile = fullfile(pname, sprintf('%s_beams.mat',prefix));
        if ~existfile(beamfile)
            fprintf('Processing beam patterns...\n')
            [beam, ref] = az_process_beams(fname1,fname2,call,array,beamfile);
        else
            fprintf('Loading precomputed beam structure...\n')
            load(beamfile)
        end

        % generate audio track
        disp(repmat('#',1,70))
        wavname = fullfile(pname, sprintf('%s.wav',prefix));
        if ~existfile(wavname)
            fprintf('Generating WAV audio track from reference data...\n')
            genAudioTrack(ref, wavname);
        else
            fprintf('WAV audio file already exists.\n')
        end

        % generate video track
        disp(repmat('#',1,70))
        for fc = plotFreq
            aviname = fullfile(pname, sprintf('%s_%dkHz.avi',prefix,round(fc*1e-3)));
            if ~existfile(aviname)
                fprintf('Generating AVI video track from beam data...\n')
                pause(5)    % mbp has trouble keeping up
                genVideoTrack(beam, ref, aviname, fc);
            else
                fprintf('AVI video track already exists.\n')
            end
        end

        % combine beam video and audio with camera video
        disp(repmat('#',1,70))
        fprinf('Combining video and audio into single AVI movie...\n')
        %TBD

        disp(repmat('#',1,70))
        fprintf('Completed processing data files with prefix "%s":\n\t%s\t%s\n',prefix,fname1,fname2)
        toc
        disp(repmat('#',1,70))
        
    catch ME
        % bail out and continue with next event
        warning('Could not process files!')
        disp(getReport(ME))
        toc
        continue
    end

    % clear memory and close figures for next run
    clear beam ref
    close all
    pause(1)
    
end

fprintf('\nCompleted processing all data files:  ')
disp(datestr(now))
diary off
