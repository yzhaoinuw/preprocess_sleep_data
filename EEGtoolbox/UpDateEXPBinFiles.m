function s=UpDateEXPBinFiles(s,BinFileName,TStart,Duration)

%s= xml2struct(XMLfilename);

if isfield(s.Animal.Acquisition.Files,'File')==0
    nFile=1;
else
    nFile=length(s.Animal.Acquisition.Files.File)+1;
end

s.Animal.Acquisition.Files.File{nFile}.FileName=sprintf('%s',BinFileName);
s.Animal.Acquisition.Files.File{nFile}.TStart=TStart;
s.Animal.Acquisition.Files.File{nFile}.Duration=Duration;

%struct2xml(s,XMLfilename);