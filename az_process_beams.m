function varargout = az_process_beams(varargin)
% AZ_PROCESS_BEAMS  analyze recorded data files and reconstruct beams
%
% az_process_beams(FNAME1,FNAME2,EVENT,ARRAY) processes the events in the
%     specified struct, EVENT, and returns the BEAM pattern data for each
%     event.  ARRAY is the struct defining array coordinates and data
%     channel mapping.
%
% az_process_beams(PREFIX,EVENT,ARRAY) uses the string PREFIX for
%     specifying the pair of SRZ file names
%
% [BEAM,SOURCE,REF,ARRAY,FREQ] = az_process_beams(...) optionally returns
%     the specified processing results
%
% [..] = az_process_beams(FNAME1,FNAME2,EVENT,ARRAY,BEAMFILE) writes beam data
%     to the specified filename, BEAMFILE
%
% Notes:
%     EVENT can be a struct from either az_detect_events or az_split_event
%     ARRAY is a struct from az_define_array
%     If not specified, BEAMFILE, will default to [PREFIX '_beams.mat']
%

warning('OFF','CALC_SPECTRUM:fs');
warning('OFF','CALC_SPECTRUM:dc');
warning('OFF','AZ_CHANINDEX:badchannel');

% plotting flags for debug
PLOT1 = 0;          % plot array channel positions
PLOT2 = 0;          % spectrogram for each raw call
PLOT3 = 0;          % 3D representation of array and source location
PLOT4 = 0;          % spectrogram for each filtered call
PLOT5 = 0;          % 3D beam surface/contour plot for each call

% optional processing modes - useful for bypassing lengthy steps if debugging
LOCMODE = false;    % localize sound source using TDOA, if true
FILTMODE = false;   % apply time-frequency filtering of each harmonic, if true
EQMODE = false;     % apply microphone calibration equalizer, if true
TLMODE = true;     % apply transmission loss correction, if true

switch nargin
    case 3
        prefix = varargin{1};
        event = varargin{2};
        array = varargin{3};
        fname1 = [prefix '_side1.srz'];
        fname2 = [prefix '_side2.srz'];
    case 4
        fname1 = varargin{1};
        fname2 = varargin{2};
        event = varargin{3};
        array = varargin{4};
    case 5
        fname1 = varargin{1};
        fname2 = varargin{2};
        event = varargin{3};
        array = varargin{4};
        beamfile = varargin{5};
        [pname,fname,ext] = fileparts(beamfile);
        if ~strcmp(ext,'.mat')
            beamfile = fullfile(pname,[fname '.mat']);
        end
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
    beamfile = [prefix 'beam.mat'];
end

% load data file, if it exists
if ischar(event)
    if exist('eventfile','file')
        fprintf('Loading events from file...')
        load(eventfile,'event');
    else
        error('Event file not found: %s',event)
    end
end


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
    source(eNum) = az_localize(array, ts, PLOT3, LOCMODE);
    
    % Realign data, separate harmonic components, and filter to remove echoes and reverb
    ts = az_filter(ts, source(eNum), array, FILTMODE);
    
    % Equalize microphone responses using calibration data
    if EQMODE
        ts = az_equalize(ts);
    else
        warning('AZ_PROCESS_BEAMS:eqmode','Bypassed microphone equilization')
    end

    % Correct data for transmission losses on each channel
    if TLMODE
        ts = az_armaloss(ts, source(eNum).rng);
    else
        warning('AZ_PROCESS_BEAMS:tlmode','Bypassed transmission loss correction')
    end
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
save(beamfile,'beam','ref','source');
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
