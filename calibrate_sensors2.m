

FORCEDET = true;

fname1 = 'caldata_side1.srz';
fname2 = 'caldata_side2.srz';

%% detect call events in file
prefix = regexp(fname1,'[_\-\ ]');      % use current filename
callfile = [fname1(1:prefix(end)) 'callmap.mat'];
if exist(callfile,'file') && ~FORCEDET
    fprintf('\nCall index file already exists!  Loading call data in "%s"...\n\n',callfile);
    load(callfile,'callmap');
else
    [callmap,hdr] = az_detect(fname1,fname2);
    fprintf('\nSaving call index to "%s"...\n\n',callfile)
    save(callfile,'callmap','hdr');
end


%% map array channels to data channels
