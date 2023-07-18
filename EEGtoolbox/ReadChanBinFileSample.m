function Data=ReadChanBinFileSample(Fs,To,NbSample,fid,CurrCh,TotnumberOfchan,Type)

nbArg=nargin;

if nbArg==6
    Type='uint16';
end

%allow to read only one channel from To to to+nbsample
ToSamp=floor(To*Fs);
% binduration=floor(Te*Fs)-floor(To*Fs);

OffsetSeek=(CurrCh-1)*2+ToSamp*TotnumberOfchan*2;

skip=(TotnumberOfchan-1)*2;

fseek(fid,OffsetSeek,'bof');
Data = fread(fid,NbSample,Type,skip);

