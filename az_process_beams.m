function varargout = az_process_beams(varargin)
% AZ_PROCESS_BEAMS  analyze Arrayzilla data files and construct beams
%
% az_process_beams(fname1,fname2) detects all calls in the files and saves
%     call map to a MAT file in the current directory
% az_process_beams(fname1,fname2,CALLS) only plots the calls in the
%     specified array
% az_process_beams(fname1,fname2,CALLS,beamfile) writes beam data to the
%     specified filename
%

warning('OFF','CALC_SPECTRUM:fs');
warning('OFF','CALC_SPECTRUM:dc');
warning('OFF','AZ_CHANINDEX:badchannel');

if ~exist('TDOA_frame','file')
    cLoc = fileparts(mfilename('fullpath'));
    addpath(fullfile(cLoc,'primary_analysis'));
end

% plotting flags
PLOT0 = 1;          % time series of detected calls
PLOT1 = 1;          % plot array channel positions
PLOT2 = 1;          % spectrogram for each raw call
PLOT3 = 0;          % 3D representation of array and source location
PLOT4 = 0;          % spectrogram for each filtered call
PLOT5 = 0;          % 3D beam surface/contour plot for each call

% force (re)detection of calls - if true, overwrites existing callmap
FORCEDET = false;

% filter mode - if true, applies time-frequency filtering around each harmonic
FILTMODE = true;

% default beam plot style
PLOTMODE = 'surf';
%PLOTMODE = 'cont';
%PLOTMODE = 'horz';
%PLOTMODE = 'vert';

CALLNUM = [];       % don't process calls unless entered, only detect them and save callmap


% prompt for filename if not entered
switch nargin
    case {0,1}
        % take first parameter as callnum if it's numeric
        if nargin == 1 && isnumeric(varargin{1})
            CALLNUM = varargin{1};
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
        
        % take specified call numbers
        CALLNUM = varargin{3};
        
    case 4
        % assign file names
        fname1 = varargin{1};
        fname2 = varargin{2};
        
        % take specified call numbers
        CALLNUM = varargin{3};
        
        % specify beam file to write
        beamfile = varargin{4};
        
    case 5
        % assign file names
        fname1 = varargin{1};
        fname2 = varargin{2};
        
        % take specified call numbers
        CALLNUM = varargin{3};
        
        % specify beam file to write
        beamfile = varargin{4};
        
        % optional plotting parameters
        FORCEDET = varargin{5};
        
    otherwise
        error('Bad number of input parameters.  Better luck next time!');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detect calls in data files and return timestamps
prefix = regexp(fname1,'[_\-\ ]');      % use current filename
callfile = [fname1(1:prefix(end)) 'callmap.mat'];
hdrfile = [fname1(1:prefix(end)) 'hdr.mat'];
if exist(callfile,'file') && ~FORCEDET
    fprintf('\nCall index file already exists!  Loading call data in "%s"...\n\n',callfile);
    load(callfile,'callmap');
else
    [callmap,hdr] = az_detect(fname1,fname2);     % detect trigger events in data header fields
    fprintf('\nSaving call index to "%s"...\n\n',callfile)
    save(callfile,'callmap');
    save(hdrfile,'hdr');
    clear hdr
end

% get call index
if CALLNUM == Inf
    callIdx = 1:length(callmap);    % assign callmap index
else
    callIdx = CALLNUM;              % otherwise use default
    if isempty(callIdx)
        return              % terminate if no calls to process
    end
end

% plot time series data on several channels
if PLOT0; plotTimeSeries(fname1, fname2, callmap(callIdx)); drawnow; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define physical array parameters
array = az_positions([19 12],[4 5]*.0254,[9/8 9/10]);

% assign channel/board mapping to array struct
array = az_channelmap(array,PLOT1);

%% preallocate structs
N = length(callIdx);
source(N).xSrc = [];
source(N).ySrc = [];
source(N).zSrc = [];
source(N).residual = [];
source(N).az = [];
source(N).el = [];
source(N).rng = [];
beam = cell(N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iterate over each call
for k = 1:length(callIdx)
    try
    
    cNum = callIdx(k);
    
    %% Convert raw recorded digital data to voltage units
    ts = az_convert(fname1,fname2,callmap(cNum),array);
    if PLOT2; plotSpecArray(array,ts); end
    
    %% Localize point sources in 3D space
    source(k) = az_localize(ts, array, PLOT3);
    
    %% Realign data, separate harmonic components, and filter to remove echoes and reverb
    ts = az_filter(ts, source(k), array, FILTMODE);
    
    %% Equalize microphone responses using calibration data
    %ts = az_equalize(ts);
    
    %% Correct data for transmission losses on each channel
    ts = az_armaloss(ts, source(k).rng);
    
    if PLOT4; plotSpecArray(array,ts); end
    
    %% Analyze frequency-content of each channel
    [fd(k),ref(k)] = az_analysis(ts,callmap(cNum));
    
    %% Interpolate beam data
    beam{k} = az_calcbeam(fd(k), array, source(k));%, 'pos', 'nearest');
    if PLOT5; plotBeamPattern(beam{k},PLOTMODE); pause; end
    
    fprintf('\n\n*** Completed processing call %d ***\n\n',cNum)
    
    %% If problem arises, issue error message and move to next call
    catch
        fprintf('ERROR: \n')
        ERR = lasterror;
        disp(ERR.message)
        
        fprintf('\n#######################################\n');
        fprintf('#####  Failed to process call %d  #####\n',cNum);
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
