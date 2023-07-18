function WriteBinFile(Fs,To,fid,NbCh,Data,Type)
%write data in bin file
%Fs sampling rate in Hz
%fid bin file handle
%NbCh number of channel
%Data Data to write nbChan x Sampl
%Type 'uint16' or 'int6'


if nargin<=5
    Type='uint16';
end

OffsetSeek=floor(To*Fs)*NbCh*2;
fseek(fid,OffsetSeek,'bof');
if strcmpi(Type,'uint16')==1
    fwrite(fid,uint16(Data),'uint16');
elseif strcmpi(Type,'int16')==1
     fwrite(fid,int16(Data),'int16');
end
