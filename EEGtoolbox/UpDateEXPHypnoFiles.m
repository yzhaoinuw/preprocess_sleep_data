function s=UpDateEXPHypnoFiles(s,HypnoFilename,TStart,Duration)

%s= xml2struct(XMLfilename);

if isfield(s.Animal.Hypnogram.Files,'File')==0
    nFile=1;
else
    nFile=length(s.Animal.Hypnogram.Files.File)+1;
end

if isstr(TStart)==0 %conversion in String
    formatOut = 'yyyy-mm-ddTHH:MM:SS.FFF';
    TStart=datestr(TStart,formatOut(1:end-4));
end

s.Animal.Hypnogram.Files.File{nFile}.FileName=sprintf('%s',HypnoFilename);
s.Animal.Hypnogram.Files.File{nFile}.TStart=TStart;
s.Animal.Hypnogram.Files.File{nFile}.Duration=Duration;

%struct2xml(s,XMLfilename);