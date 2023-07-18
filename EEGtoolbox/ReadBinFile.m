function [Data]=ReadBinFile(Fs,To,Te,fid,NbCh,skip,Type)

nbArg=nargin;

if nbArg==6
    Type='uint16';
end

%OffsetSeek=(Nch-1)*2+To*Fs*NbCh*2;
OffsetSeek=floor(To*Fs)*NbCh*2;
binduration=floor(Te*Fs)-floor(To*Fs);
fseek(fid,OffsetSeek,'bof');


Format=sprintf('%d*%s',NbCh,Type);

if skip>1
    skipValue=NbCh*2*(skip-1);
    
    Data = fread(fid,[NbCh,binduration/skip],Format,skipValue);
else
    Data = fread(fid,[NbCh,binduration],Format);
end

