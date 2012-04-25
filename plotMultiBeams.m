function plotMultiBeams(B,varargin)
% PLOTMULTIBEAMS  plots multiple beam patterns over time
%
% PLOTMULTIBEAMS(B) iterates over each beam in the structure, B, to show
% the progression over time.  User can select the parameters to iterate over
% between azimuth, elevation, or frequency.
%
% Beam is an Mx1 struct containing the 3-dimensional array of complex beam
% data, and indices for each dimension.
%
% beam.Z - i x j x k complex matrix azimuth, elevation, and frequency
% beam.y - i-element vector of elevation indices
% beam.x - j-element vector of azimuth indices
% beam.f - k-element vector of frequency indices
% beam.Y - i x j matrix of elevation for each element in Z
% beam.X - i x j matrix of azimuth for each element in Z
% 
% To generate X and Y use:
% [beam.X, beam.Y] = meshgrid(beam.x, beam.y);


fprintf('\n\n*****************************************\n')
fprintf('Plotting beam patterns\n\n')

PLOTMODE = 'freq'; %'az','el'
if nargin > 1
    PLOTMODE = varargin{1};
end

