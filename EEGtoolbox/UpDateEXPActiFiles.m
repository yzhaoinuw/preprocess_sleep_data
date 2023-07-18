function s=UpDateEXPActiFiles(s,Acti,CurrVideo)

%creation of the field video if doesn't exist
if isfield(s.Animal.Videos.Video{CurrVideo},'Actimetry')==0
    nActi=1;
   
    s.Animal.Videos.Video{CurrVideo}.Actimetry.Name=sprintf('%s',Acti.Name);
    s.Animal.Videos.Video{CurrVideo}.Actimetry.SamplingRate=s.Animal.Videos.Video{CurrVideo}.SamplingRate;
    s.Animal.Videos.Video{CurrVideo}.Actimetry.Enable='true';
    s.Animal.Videos.Video{CurrVideo}.Actimetry.Files=[];

    
end


%check how many files tehre is the field Files
if isfield(s.Animal.Videos.Video{CurrVideo}.Actimetry,'Files')==1 && isfield(s.Animal.Videos.Video{CurrVideo}.Actimetry.Files,'File')==1
   nbfile=length(s.Animal.Videos.Video{CurrVideo}.Actimetry.Files.File);
else
   nbfile=0;
end
s.Animal.Videos.Video{CurrVideo}.Actimetry.Files.File{nbfile+1}.FileName=sprintf('%s',Acti.FileName);
s.Animal.Videos.Video{CurrVideo}.Actimetry.Files.File{nbfile+1}.TStart=Acti.TStart;
s.Animal.Videos.Video{CurrVideo}.Actimetry.Files.File{nbfile+1}.Duration=round(Acti.Duration*1000);


