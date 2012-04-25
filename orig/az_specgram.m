close all

% convert SRZ data to MAT format for analysis
%pname = 'data/20110623';
%fname = 'piston_side1.srz';   %'chris_hand.srz';  % 'jim_clapping.srz'; %'sideA_1second.srz';%
%fname = 'apollo_side1.srz';  %'torpi_side1.srz';

%fname = 'piston_side1.srz';
%idx = 6.9e4:7.05e4;  %20110629/piston_side1
%fname = 'piston_side2.srz';
%idx = 6.0e4:6.15e4;

%fname = '20110715_test_bats_side_1.srz';
%idx = 5.34e5:5.36e5;
%idx = 5.08750e5:5.10250e5;
%ch=;
%fname = '20110715_test_bats_side_2.srz';
%idx = idx-300;
% idx = 5.34e5:5.36e5;
% idx = 4.0e4:4.152e4;

% fname = 'torpi_side1.srz';
% idx = 8.825e4:8.96e4;
%fname = 'torpi_side2.srz';
%idx = 8.0e4:8.135e4;

%fname = '20110715_test2_side1.srz';
%idx = 1:1e6;
%idx = 4.32e5:4.36e5;
%fname = '20110715_test2_side2.srz';
%ch = 67;

%fname = '20110718_speaker2_side1.srz';
%fname = '20110718_speaker2_side2.srz';

%fname = '20110718_speaker_angled_side1.srz';
%fname = '20110718_speaker_angled_side2.srz';
%fname = '20110718_speaker_newmics_side1.srz';
fname = '20110718_half_side1.srz';

idx = 1:236660;
%idx = 700:2000;
%idx = (700:2500) + 236660.7*.02;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIDE 1
% quad 1
% ch = [1 29 57 85 31 59 87 4 61 89 6 34 91 8 36 64];

% quad 2
%ch = [2 30 58 86 32 60 88 5 62 90 7 35 92 9 37 65];

% quad 3
%ch = [10 38 66 94 40 68 96 13 70 98 15 43 17 45 73 101];

% quad 4
%ch = [11 39 67 95 41 69 97 14 71 99 16 44 18 46 74 102];

% quad 5
%ch = [75 103 20 48 22 50 78 106 80 108 25 53 110 27 55 83];

% quad 6
%ch = [76 104 21 49 23 51 79 107 81 109 26 54 111 28 56 84];

% lcol
ch = [3 33 63 93 12 42 100 47 105 52 72 19 77 24 82 112];

%%%%%%%%%%%%%%%%%%%%%%%%
% SIDE 2
% quad 1
% ch = [111 83 55 27 52 24 109 81 106 78 50 22 47 19 104 76];

% quad 2
%ch = [112 84 56 28 53 25 110 82 107 79 51 23 48 20 105 77];

% quad 3
%ch = [101 73 45 17 42 14 99 71 96 68 40 12 9 94 66 38];

% quad 4
ch = [102 74 46 18 43 15 100 72 97 69 41 13 10 95 67 39];

% quad 5
%ch = [35 7 92 64 61 33 5 90 87 59 31 3 86 58 30 2];

% quad 6
% ch = [36 8 93 65 62 34 6 91 88 60 32 4 85 57 29 1];

% rcol
%ch = [54 108 49 103 44 98 11 37 63 89 26 80 21 75 16 70];

%ch = [10 95 67 36 8 93 62 34 6];  % lower grid
%ch = [26 80 21 75 16 70];
%idx = 1:300e3;
%idx = 127e3:128e3;
%idx = 117e3:118e3;
%ch = [80 25 110 104 101 18 14 71 43 97 69 38 39 33 62 4];  %bad mic channels
%ch = [112 82 79 108 13 35 92 37];  % bad fai channels

[ts.data, hdr] = read_SRZ(fname,idx,ch);        % read in data set

ts.data = ts.data * 5*2^-16;        % convert to Volts
ts.data = ts.data - 2.5;            % remove DC offset

figure(999)
plot(ts.data)
grid on

%%% analysis

% load board A
for i = 1:length(ch)
    [S,F,T,P] = spectrogram(ts.data(:,i),256,240,256,hdr.fs,'yaxis');
    
    figure(ch(i))
    set(gcf,'MenuBar','none')
    set(gcf,'ToolBar','figure')
    imagesc(T,F,10*log10(abs(P)));
    set(gca,'YDir','normal')
    set(gca,'CLim',[-120 -50]);
    colorbar
end

tilefigs(ch,4,4);

% load board B

