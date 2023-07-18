function Data=ReadChanBinFile(Fs,To,Te,fid,CurrCh,TotnumberOfchan,Type)
%read only one channel from To to Te
%Fs: sampling rate
%T0 in sec
%Te in sec
%fid handles of the bin file to read
%CurrCh Chan to read
%TotnumberOfchan number of channel in the binary file
%Type 'uint16' or 'int16' deafaul uint16 could be ignored

nbArg=nargin;

if nbArg==6
    Type='uint16';
end

ToSamp=floor(To*Fs);
binduration=floor(Te*Fs)-floor(To*Fs);

OffsetSeek=(CurrCh-1)*2+ToSamp*TotnumberOfchan*2;

skip=(TotnumberOfchan-1)*2;

fseek(fid,OffsetSeek,'bof');
Data = fread(fid,binduration,Type,skip);

