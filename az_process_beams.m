function varargout = az_process_beams(varargin)
% AZ_PROCESS_BEAMS  analyze Arrayzilla data files and construct beams
%
% BEAM = az_process_beams(FNAME1,FNAME2,EVENT) processes the events in the
%     specified struct and returns the BEAM pattern data for each event
%
% [BEAM,SOURCE,REF,ARRAY,FREQ] = az_process_beams(...) optionally returns
%     the specified additional processing results
%
% [..] = az_process_beams(FNAME1,FNAME2,EVENTS,BEAMFILE) writes beam data
%     to the specified filename, BEAMFILE
%

warning('OFF','CALC_SPECTRUM:fs');
warning('OFF','CALC_SPECTRUM:dc');
warning('OFF','AZ_CHANINDEX:badchannel');

% if ~exist('TDOA_frame','file')
%     cLoc = fileparts(mfilename('fullpath'));
%     addpath(fullfile(cLoc,'primary_analysis'));
% end

% plotting flags
PLOT1 = 0;          % plot array channel positions
PLOT2 = 1;          % spectrogram for each raw call
PLOT3 = 0;          % 3D representation of array and source location
PLOT4 = 1;          % spectrogram for each filtered call
PLOT5 = 0;          % 3D beam surface/contour plot for each call

% filter mode - if true, applies time-frequency filtering around each harmonic
FILTMODE = false;

% prompt for filename if not entered
switch nargin
    case {0,1}
        % take first parameter as eventnum if it's numeric
        if nargin == 1 && isnumeric(varargin{1})
            EVENTNUM = varargin{1};
        % take first parameter as FORCEDET if boolean
        elseif nargin == 1 && islogical(varargin{1})
            FORCEDET = varargin{1};
        end
        
        % prompt user for data file
        [fname, pname] = uigetfile( { ...
            '*.srz','Recorder Data (*.srz)'; ...
            '*.*',  'All Files (*.*)'}, ...
            'Select Side 1 data file', ...
            'MultiSelect', 'on');
        
        if iscell(fname) && length(fname) == 2
            % assign 2nd filename if multiple selected
            fname1 = fullfile(pname,fname{1});
            fname2 = fullfile(pname,fname{2});
        else
            % otherwise get 2nd file name
            fname1 = fullfile(pname,fname);
            [fname, pname] = uigetfile( { ...
                '*.srz','Recorder Data (*.srz)'; ...
                '*.*',  'All Files (*.*)'}, ...
                'Select Side 1 data file');
            fname2 = fullfile(pname,fname);
        end
        clear fname pname
        
    case 2
        % assign filenames
        fname1 = varargin{1};
        fname2 = varargin{2};
        
    case 3
        % assign file names
        fname1 = varargin{1};
        fname2 = varargin{2};
        
        % specify beam file to write
        beamfile = varargin{3};
        
    otherwise
        error('Bad number of input parameters.  Better luck next time!');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detect trigger events in data files and return index
prefix = regexp(fname1,'[_\-\ ]');      % use current filename
prefix = fname1(1:prefix(end));
eventfile = [prefix 'events.mat'];
hdrfile = [prefix 'hdr.mat'];
if exist(eventfile,'file') && ~FORCEDET
    fprintf('Loading events from file...')
    load(eventfile,'events');                               % load data file, if exists
    fprintf(' Done\n')
else
    events = az_detect_events(fname1,fname2,eventfile,hdrfile);    % detect trigger events in data header fields
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define physical array parameters
array = az_positions([19 12],[4 5]*.0254,[9/8 9/10]);

% assign channel/board mapping to array struct
array = az_channelmap(array,PLOT1);     % MOVE BAD CHANNEL SPECIFICATION, DETECTION & REDUCTION TO LATER

%% preallocate structs
N = length(eIdx);

[source(1:N).xSrc] = deal([]);
[source(1:N).ySrc] = deal([]);
[source(1:N).zSrc] = deal([]);
[source(1:N).residual] = deal([]);
[source(1:N).az] = deal([]);
[source(1:N).el] = deal([]);
[source(1:N).rng] = deal([]);

beam = cell(N,1);

[ref(1:N).eNum] = deal([]);
[ref(1:N).t0] = deal([]);
[ref(1:N).t1] = deal([]);
[ref(1:N).tlen] = deal([]);
[ref(1:N).ch] = deal([]);
[ref(1:N).data] = deal([]);
[ref(1:N).fs] = deal([]);
[ref(1:N).done] = deal(false);
[ref(1:N).error] = deal([]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iterate over each event
%cNum = 0;      % initialize number of calls detected in all events
for eNum = eIdx
    try

%    % look for multiple calls in each event; if found, iterate over each one
%    callfile = sprintf('%scalls_%d',prefix, eNum);
%    calls = az_split_event(fname1,fname2,events(eNum),array,callfile);
    
%     % iterate over each call in event
%     for m = 1:numel(calls)
%         cNum = cNum + 1;
%         
%         if numel(calls) > 1
%             fprintf('\n********************************\n')
%             fprintf('*** Processing call %d of %d ***\n\n',m,numel(calls))
%         end

        %% Convert raw recorded digital data to voltage units
        ts = az_convert(fname1,fname2,calls(m),array);
        if PLOT2; plotSpecArray(array,ts); end

        %% Localize point sources in 3D space
        source(cNum) = az_localize(ts, array, PLOT3);

        %% Realign data, separate harmonic components, and filter to remove echoes and reverb
        ts = az_filter(ts, source(cNum), array, FILTMODE);

        %% Equalize microphone responses using calibration data
        %ts = az_equalize(ts);

        %% Correct data for transmission losses on each channel
        ts = az_armaloss(ts, source(cNum).rng);
        if PLOT4; plotSpecArray(array,ts); end

        %% estimate bulk parameters
        ref{cNum} = az_estimate_params(ts,calls(m));

        %% Analyze frequency-content of each channel
        fd(cNum) = az_analysis(ts);

        %% Interpolate beam data
        beam{cNum} = az_calcbeam(fd(cNum), array, source(cNum));%, 'pos', 'nearest');

%        if PLOT5; plotBeamPattern(beam{cNum},60e3,PLOTMODE); pause; end
%    end
    
    fprintf('\n***************************************\n')
    fprintf('*** Completed processing event %.3d ***\n',eNum)
    fprintf('***************************************\n\n')
    ref{cNum}.done = true;             % set done flag

    %% If problem arises, issue error message and move to next event
    catch ME
        ref{cNum}.error = ME;
        disp(getReport(ME));
        fprintf('#######################################\n');
        fprintf('#####  Failed to process event %d  #####\n',eNum);
        fprintf('#######################################\n\n');
    end
end


% save data
if ~exist('beamfile','var')
    beamfile = [fname1(1:prefix(end)) 'beams.mat'];
end

% write data to new file, if already exists
if exist(beamfile,'file')
    k = 1;
    beamfile = sprintf('%s_%.2d.mat',beamfile(1:end-4),k);
    while exist(beamfile,'file')
        k = k+1;
        beamfile = sprintf('%s_%.2d.mat',beamfile(1:end-7),k);
    end
end

% save to beamfile
save(beamfile,'beam','ref','source','array');
fprintf('Saving beam data to %s\n',beamfile);


%% send appropriate output
switch nargout
    case 0
    case 1
        varargout{1} = beam;
    case 2
        varargout{1} = beam;
        varargout{2} = ref;
    case 3
        varargout{1} = beam;
        varargout{2} = ref;
        varargout{3} = source;
    case 4
        varargout{1} = beam;
        varargout{2} = ref;
        varargout{3} = source;
        varargout{4} = array;
    case 5
        varargout{1} = beam;
        varargout{2} = ref;
        varargout{3} = source;
        varargout{4} = array;
        varargout{5} = fd;
end
