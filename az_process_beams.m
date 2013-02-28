function varargout = az_process_beams(varargin)
% AZ_PROCESS_BEAMS  analyze Arrayzilla data files and construct beams
%
% az_process_beams(FNAME1,FNAME2,EVENT) processes the events in the
%     specified struct and returns the BEAM pattern data for each event
%
% az_process_beams(PREFIX,EVENT) uses the string PREFIX for SRZ file names
%
% [BEAM,SOURCE,REF,ARRAY,FREQ] = az_process_beams(...) optionally returns
%     the specified processing results
%
% [..] = az_process_beams(FNAME1,FNAME2,EVENT,BEAMFILE) writes beam data
%     to the specified filename, BEAMFILE
%
% EVENT can be a struct from either az_detect_events or az_split_event

warning('OFF','CALC_SPECTRUM:fs');
warning('OFF','CALC_SPECTRUM:dc');
warning('OFF','AZ_CHANINDEX:badchannel');

% plotting flags
PLOT1 = 0;          % plot array channel positions
PLOT2 = 0;          % spectrogram for each raw call
PLOT3 = 0;          % 3D representation of array and source location
PLOT4 = 0;          % spectrogram for each filtered call
PLOT5 = 0;          % 3D beam surface/contour plot for each call

% filter mode - if true, applies time-frequency filtering around each harmonic
FILTMODE = false;

switch nargin
    case 2
        prefix = varargin{1};
        event = varargin{2};
        
        fname1 = [prefix '_side1.srz'];
        fname2 = [prefix '_side2.srz'];
    case 3
        fname1 = varargin{1};
        fname2 = varargin{2};
        event = varargin{3};
    case 4
        fname1 = varargin{1};
        fname2 = varargin{2};
        event = varargin{3};
        
        % specify beam file to write
        beamfile = varargin{4};
        
    otherwise
        error('Bad number of input parameters.  Better luck next time!');
end

% get filename prefix from current filename
if ~exist('prefix','var')
    prefix = regexp(fname1,'[_\-\ ]');
    prefix = fname1(1:prefix(end));
end

% assign beamfile name if not specified
if ~exist('beamfile','var')
    beamfile = [prefix 'beam'];
end

% load data file, if it exists
if ischar(event)
    if exist(eventfile,'file')
        fprintf('Loading events from file...')
        load(eventfile,'event');
    else
        error('Event file not found: %s',event)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define physical array parameters
array = az_positions([19 12],[4 5]*.0254,[9/8 9/10]);

% assign channel/board mapping to array struct
array = az_channelmap(array,PLOT1);     % MOVE BAD CHANNEL SPECIFICATION, DETECTION & REDUCTION TO LATER

%% preallocate structs
N = numel(event);

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
for eNum = 1:N
    try

    % Convert raw recorded digital data to voltage units
    ts = az_convert(fname1,fname2,event(eNum),array);
    if PLOT2; plotSpecArray(array,ts); end

    % Localize point sources in 3D space
    source(eNum) = az_localize(ts, array, PLOT3);

    % Realign data, separate harmonic components, and filter to remove echoes and reverb
    ts = az_filter(ts, source(eNum), array, FILTMODE);

    % Equalize microphone responses using calibration data
    %ts = az_equalize(ts);

    % Correct data for transmission losses on each channel
    ts = az_armaloss(ts, source(eNum).rng);
    if PLOT4; plotSpecArray(array,ts); end

    % estimate bulk parameters
    ref(eNum) = az_estimate_params(ts,event(eNum));

    % Analyze frequency-content of each channel
    fd(eNum) = az_analysis(ts);

    % Interpolate beam data
    beam{eNum} = az_calcbeam(fd(eNum), array, source(eNum));%, 'pos', 'nearest');
    if PLOT5; plotBeamPattern(beam{eNum},60e3,PLOTMODE); pause; end
    
    fprintf('\n\n**************************************\n')
    fprintf('*** Completed processing event %.3d ***\n',eNum)
    fprintf('**************************************\n\n')
    ref(eNum).done = true;             % set done flag

    % If problem arises, issue error message and move to next event
    catch ME
        ref(eNum).error = ME;
        disp(getReport(ME));
        fprintf('#######################################\n');
        fprintf('#####  Failed to process event %d  #####\n',eNum);
        fprintf('#######################################\n\n');
    end
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
