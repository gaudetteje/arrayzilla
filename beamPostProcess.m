% compute and plot bulk parameters from beam data

fname = 'dori2_worm_beams_1-647.mat';

prefix = 'Dori2';



% main response axis location (x vs. y)
zpk(f) = max(max(Z));
[xpk(f), ypk(f)] = find(Z==zpk(f));

% pulse length vs. time
figure; plot([ref.t0],1e3*[ref.tlen],'.')
xlabel('Time (seconds)')
ylabel('Pulse length (ms)')
title(sprintf('Pulse Duration - %s',prefix))
% save to matlab fig
% save to eps figure