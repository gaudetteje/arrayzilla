function array = az_define_array(varargin)
% AZ_DEFINE_ARRAY  defines array coordinates and board/channel mapping
% that corresponds to data channels in the SRZ raw recorder files
%
% ARRAY = az_define_array() returns a struct, ARRAY, containing coordinates
%     'xPos' and 'yPos' along with the board and channel mappings, 'bd'
%     and 'ch'.  Known bad channels are explicitly excluded from the
%     struct, but listed in the fields 'badCh1' and 'badCh2'.
%
% ARRAY = az_define_array(FNAME) saves the array struct definition to the
%     file, FNAME
%
% ARRAY = az_define_array(..,true) plots the 2D microphone positions in
%     array coordinates with the Z axis corresponding to channel number
%     for identification
%
% Note:  Channel and board numbers are hard coded values and must be
%        updated manually if array configuration is ever changed.


arrayfile = [];
PLOTFLAG = false;
switch nargin
    case 0
    case 1
        if ischar(varargin{1})
            arrayfile = varargin{1};
        else
            PLOTFLAG = varargin{1};
        end
    case 2
        arrayfile = varargin{1};
        PLOTFLAG = varargin{2};
end

% first, define array spacing and grid
array = az_positions([19 12],[4 5]*.0254,[9/8 9/10]);

% bump in first column of sensors to match hardware
array.xPos(array.xPos < 0) = 0.0127;

% Assign each coordinate to channel number in recorded data (use 0 for no connection)
array.bd = [0 ones(1,9) 2*ones(1,8) 0 ...
        ones(1,10) 2*ones(1,9) ...
        ones(1,10) 2*ones(1,9) ...
        ones(1,10) 2*ones(1,9) ...
        ones(1,10) 2*ones(1,9) ...
        ones(1,10) 2*ones(1,9) ...
        ones(1,9) 2*ones(1,10) ...
        ones(1,9) 2*ones(1,10) ...
        ones(1,9) 2*ones(1,10) ...
        ones(1,9) 2*ones(1,10) ...
        ones(1,9) 2*ones(1,10) ...
        0 ones(1,8) 2*ones(1,9) 0];
%surf(array.xPos(1:19), array.yPos(1:19:end), reshape(array.bd,19,12)');   % for sanity check

array.ch = [  0 110  27  55  83 111  28  56  84 112  86  58  30   2  85  57  29   1   0 ...
         52  80 108  25  53  81 109  26  54  82  87  59  31   3  88  60  32   4  89 ...
        105  22  50  78 106  23  51  79 107  24  61  33   5  90  62  34   6  91  63 ...
         47  75 103  20  48  76 104  21  49  77  35   7  92  64  36   8  93  65  37 ...
        100  17  45  73 101  18  46  74 102  19   9  94  66  38  10  95  67  39  11 ...
         42  70  98  15  43  71  99  16  44  72  96  68  40  12  97  69  41  13  98 ...
         12  40  68  96  13  41  69  97  14  70  42  14  99  71  43  15 100  72  44 ...
         93  10  38  66  94  11  39  67  95  16 101  73  45  17 102  74  46  18 103 ...
         63  91   8  36  64  92   9  37  65  75  47  19 104  76  48  20 105  77  49 ...
         33  61  89   6  34  62  90   7  35  21 106  78  50  22 107  79  51  23 108 ...
          3  31  59  87   4  32  60  88   5  80  52  24 109  81  53  25 110  82  54 ...
          0   1  29  57  85   2  30  58  86  26 111  83  55  27 112  84  56  28   0 ...
        ];
%reshape(array.ch,19,12)';      % for sanity check

% placeholders for corner elements
idx = [1 19 210 228];

% known bad on FAI
badCh1 = [109 73]; % 74 69 50
badCh2 = [102];%2 8 9 13 35 40 102 79 108 82 112 36 37];

% known bad on microphones?
badCh1 = [badCh1 9 80 111 58 105];% 27 29 55 58 80 83 105 110 111]; %18 44 92 95];
badCh2 = [badCh2 25 53 17 71 84 112 80]; %1 2 8 18 25 28 29 53 56 57 80 88 84 85];% 20 18 93 80 13 107 25 53 28 77];% 82 23 108 40 41 13 98]; % 88 93 37];%2 8 9 13 36 40 82 96 102];
badCh1 = [];
badCh2 = [];

% locate matching channel numbers
idx = [idx az_chanindex([badCh1 badCh2], [ones(size(badCh1)) 2*ones(size(badCh2))], array)];

% remove offending channels from array structure
array.xPos(idx) = [];
array.yPos(idx) = [];
array.bd(idx) = [];
array.ch(idx) = [];

% save channel numbers that were removed
array.badCh1 = badCh1;
array.badCh2 = badCh2;

% save results to file
if ~isempty(arrayfile)
    fprintf('Saving array definition to "%s"...\n',arrayfile)
    save(arrayfile,'array');
end

% plot results
if PLOTFLAG
    figure
    plot3(array.xPos,array.yPos,array.ch,'r.');
    title('Microphone positions for Arrayzilla')
    xlabel('X position (m)')
    ylabel('Y position (m)')
    grid on
    hold on
    view(0,90)

    plot(0,0,'k+')            % origin of absolute array
    axis equal
    dx = max(0, max(diff(sort(array.xPos))));
    dy = max(0, max(diff(sort(array.yPos))));
    axis([min(array.xPos - dx/2) max(array.xPos + dx/2) ...
          min(array.yPos - dy/2) max(array.yPos + dy/2)])
    
    drawnow
end
