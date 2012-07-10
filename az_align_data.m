function az_align_data(varargin)
% AZ_ALIGN_DATA  verifies and corrects packet alignment in raw binary data
%
% az_align_data(fname) scans through data looking for "Form" at every 256
%     bytes.  Properly formatted data contains this header marker.  When
%     absent, user is prompted for what action to take.
% az_align_data(fname,offset) starts at the specified offset in bytes
% az_align_data() prompts user for a file
%
% Occasional data corruption may occur, which is unfortunately normal for
% the flashed firmware version.  This function will remove any malformed
% data and realign the raw data packets every 256 bytes.  Care should be
% taken to backup your data and verify proper operation before trusting
% this utility!!!


% define default number of 256 byte packets per block
blockSize = 131072; %2048;

offset = 0;


% parse input parameters
switch nargin
    case 0
        [fname, pname] = uigetfile( { ...
            '*.srz','Recorder Data (*.srz)'; ...
            '*.*',  'All Files (*.*)'}, ...
            'Select data file');
        fname = fullfile(pname,fname);
        
    case 1
        fname = varargin{1};
    case 2
        fname = varargin{1};
        offset = varargin{2};
        
    otherwise
        error('Incorrect number of input parameters')
end

% verify file exists
if ~exist(fname,'file')
    error('SYNC:FNF','File not found!  "%s"',fname)
end

% get file statistics
fstats = dir(fname);
nPackets = fstats.bytes/256;
fprintf('File found containing %g MB (%.1f data packets)\n',fstats.bytes/1024^2,nPackets)

% open file for rw
fid = fopen(fname,'r+');

% start at offset
if offset
    res = fseek(fid,offset,-1);
    if res == 0
        fprintf('Skipping to 0x%x (%.2f MB)\n',offset,offset/1024^2)
    else
        error('Can''t seek to specified location');
    end
end




tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iterate over blocks until end of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
code = uint32(hex2dec('6d726f46')); % compare "Form" with uint32 result (force into little endian)
blockNum = 0;   % init
fprintf('Processing %u blocks of size %d KB\n',ceil(nPackets/blockSize),blockSize*256/1024)
%fprintf('  0');
while ~feof(fid)
    
    res = fread(fid,blockSize,'*uint32',252);  % read "Form" every 256 bytes
    
    if any(res ~= code)
        fprintf('\n');
        
        % locate start of last good packet (or start of current block)
        idx = offset + 256*blockSize*blockNum + max(256*(find(res~=code,1)-2), 0);
        warning('SYNC:misalign','Misalignment detected in "%s"',fname)
        
        % prompt user for action
        fprintf('\nMisaligned packet after location 0x%X (%g%%) - Packet %u of %u\n', ...
            idx, 100*idx/fstats.bytes, floor(idx/256), floor(fstats.bytes/256))
        key = input('How should we fix the problem?\n  [I]gnore, [a]bort, [f]ix, or [d]elete all remaining samples:  ','s');
        if isempty(key); key = 'i'; end
        
        % modify file contents
        switch(lower(key(1)))
            % delete remaining bytes in file
            case 'd'
                key = input(sprintf('\nThis will remove the last %u bytes remaining (from 0x%x) in the file!\nAre you absolutely sure? [y/N]:  ',fstats.bytes-idx,idx),'s');
                if isempty(key); key = 'n'; end
                if strcmpi(key(1),'y')
                    java.io.RandomAccessFile(fname,'rw').setLength(idx);
                    break
                else
                    fprintf('Doing nothing...\n\n')
                    break
                end
                
            % modify file to concatenate good packets
            case 'f'
                
                % go back and search every byte for next header
                if fseek(fid,idx,-1) == 0
                    while ~feof(fid)
                        res = fread(fid,1024,'*uint8');
                        idx2 = strfind(res','Form');

                        % when found, start looking for next data misalignment
                        if numel(idx2) > 1

                            % mark location of last and next packet header
                            dstloc = idx;
                            srcloc = ftell(fid) - numel(res) + idx2(2) - 1;
                            fprintf('Found valid data at 0x%x\n',srcloc)
                            
                            % look for valid CRC checksums, otherwise keep looking
%                             for n = 2:numel(idx2);
%                                 if ~check_SRZ_CRC32(fname1,srcloc+252)
%                                     continue
%                                 else
%                                     srcloc = srcloc + idx2(n);  % need to prevent overflow
%                                     break
%                                 end
%                             end
                            
                            % start looking recursively at this offset (fix data backwards to avoid redundant writes)
                            fprintf('*******************\n');
                            az_align_data(fname,srcloc);
                            fprintf('*******************\n');
                            
                            fprintf('Deleting data between 0x%x and 0x%x (%u bytes)\n',idx,srcloc,srcloc-idx-1)
                            fprintf('Moving %g MB in blocks of size %d KB\n', (fstats.bytes-srcloc)/1024^2, blockSize*256/1024)
                            
                            % correct the error when nested call returns
                            fseek(fid,srcloc,-1);
                            while 1
                                % read data block starting at 'nextloc'
                                res = fread(fid,blockSize*64,'*uint32');
                                srcloc = ftell(fid);
                                
                                % write data block to 'currloc'
                                fseek(fid,dstloc,-1);
                                cnt = fwrite(fid,res,'uint32');
                                dstloc = ftell(fid);
                                fseek(fid,srcloc,-1);
                                
                                % verify successful write
                                assert(cnt == numel(res),sprintf('Read %u bytes, but only wrote %d!',numel(res),cnt))
                                
                                fprintf('.');

                                % reach end of file?
                                if feof(fid) || (fstats.bytes == srcloc) || (cnt < 64)
                                    break
                                end
                                
                            end
                            
                            % resize file to remove unused bytes
                            java.io.RandomAccessFile(fname,'rw').setLength(dstloc);
                            
                            fprintf('  Done!\n')
                            fprintf('... Deleted %d bytes after 0x%x ...\n',fstats.bytes-dstloc,dstloc)
                            break
                        end
                    end
                    
                    % if 'Form' never found, just delete data
                    if ~exist('dstloc','var')
                        fprintf('No valid data found after 0x%x\n',idx)
                        key = input(sprintf('Delete the last %u bytes remaining in the file? [y/N]:  ',fstats.bytes-idx),'s');
                        if isempty(key); key = 'n'; end
                        if strcmpi(key(1),'y')
                            java.io.RandomAccessFile(fname,'rw').setLength(idx);
                        else
                            fprintf('Doing nothing...\n\n')
                        end
                    end
                    break       % stop checking remaining blocks (remaining bytes were checked recursively)

                else
                    error('Couldn''t rewind file!')
                end

                % overwrite from last valid packet

            case 'a'
                fprintf('Aborting...\n\n');
                break

            otherwise
                fprintf('Doing nothing...\n\n')
        end
    end
    
    blockNum = blockNum+1;      % increment block count
    fprintf('.');
    %fprintf('%s%3.d',char([8 8 8]),blockNum);
    if ~rem(blockNum,80), fprintf('\n'); end    % wrap text on command line
end
fprintf('\n')
toc
fclose(fid);

fprintf('\nDone processing "%s"!\n\n',fname)
