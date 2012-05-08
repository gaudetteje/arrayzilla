function [idx,bad] = az_chanindex(ch,bd,array)
% AZ_CHANINDEX  returns the data structure index where the desired channels
% exist
%
% idx = AZ_CHANINDEX(ch,bd,array) takes the array struct and vectors ch and
%       bd and returns the column index where those channels are located.
% [idx,bad] = AZ_CHANINDEX(...) also returns an index of unlocated
%       channels, if any were found
%
% There is a one-to-one mapping for each (ch,bd) pair.

assert(numel(ch) == numel(bd), 'ch size must match bd size')

[chMtx, chIdx] = sortrows([array.bd; array.ch]');

idx = zeros(1,length(ch));
bad = [];
for n = 1:length(ch)
     arrayIdx = chIdx(chMtx(:,2) == ch(n) & chMtx(:,1) == bd(n));
     if arrayIdx
         idx(n) = arrayIdx;
     else
         warning('AZ_CHANINDEX:badchannel',sprintf('Channel %d on board %d does not exist in array!  Index will be zero.',ch(n),bd(n)))
         bad = [bad n];
     end
end

