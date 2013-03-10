%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  This file gives an overview of the arrayzilla post-processing steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% processing and data analysis functions

% az_align_data(fname1) - validates and fixes data misalignment in 1st file
% az_align_data(fname2) - validates and fixes data misalignment in 2nd file
%
% az_detect_events(fname1,fname2) - detects discontinuous recording events
% using header field information.  Discontinuous events may be due to
% explicitly triggered data, power cycling the digital boards on the array,
% or data corruption and/or hardware reset.
%
% az_split_events(fname1,fname2,events) - detects echolocation calls in
% each event by searching through time series data


%%% plotting and visualization functions

% plotArrayPoints(array) - plots geometry of array elements and
% localized source position.  Angular coordinates are mapped onto a sphere
% around the source position.
%
% plotTimeSeries(fname1,fname2,events) - plots the time series data from
% several representative channels.  events may be a single event or call
% structure or multiple events or call structures.
%
% plotAuxChan(hdr,events) - plots the time series data on the auxiliary
% channel for all events.  events may be a single event or call structure
% or multiple events or call structures.
%
% plotBeamPattern(beam{idx}) - plots the processed beam pattern with a
% variety of options (see help file)