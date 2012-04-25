function az_processtrial(fname,varargin)
% main analysis function for arrayzilla trials. 
%
%


TRIGGERMIC=1;


% if no source position is given: 
%   find calls on single specified channel
%   localize calls
% end
%
% beamform on source position (track) and find calls on beamformed track
% filter for attenuation and calcualte beampattern
% perform further analysis
%