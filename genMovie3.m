%function az_genmovie3
% AZ_GENAVI  compiles beam pattern images and audio data into an AVI file


% default parameters
D = 5;                         % slow audio playback by factor of D
tBuf = 0.1;                     % add short buffer to beginning
nbits = 16;
prefix = 'starbuck_test';

fs = 236660;

