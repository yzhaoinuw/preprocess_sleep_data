function [s,nVideo]=UpDateEXPVideoFiles(s,Video)

%creation of the field video if doesn't exist
if isfield(s.Animal.Videos,'Video')==0
    nVideo=1;
    s.Animal.Videos.Video{nVideo}.Name=sprintf('%s',Video.Name);
    s.Animal.Videos.Video{nVideo}.SamplingRate=Video.FrameRate;
    s.Animal.Videos.Video{nVideo}.Enable='true';
    s.Animal.Videos.Video{nVideo}.Files=[];
    
  
else %if more than one video
    %check if this videoname already exist
    for n=1:length(s.Animal.Videos.Video)
     listvideoname{n}=s.Animal.Videos.Video{n}.Name;
    end
     
     nVideo=find(strcmp(listvideoname,Video.Name)==1);
     if isempty(nVideo) %if video does not exist
        nVideo= length(s.Animal.Videos.Video)+1;
        s.Animal.Videos.Video{nVideo}.Name=sprintf('%s',Video.Name);
        s.Animal.Videos.Video{nVideo}.SamplingRate=Video.FrameRate;
        s.Animal.Videos.Video{nVideo}.Enable='true';
        s.Animal.Videos.Video{nVideo}.Files=[];
     end

    
end


%check how many files tehre is the field Files
if isfield(s.Animal.Videos.Video{nVideo},'Files')==1 && isfield(s.Animal.Videos.Video{nVideo}.Files,'File')==1
   nbfile=length(s.Animal.Videos.Video{nVideo}.Files.File);
else
   nbfile=0;
end
s.Animal.Videos.Video{nVideo}.Files.File{nbfile+1}.FileName=sprintf('%s',Video.FileName);
s.Animal.Videos.Video{nVideo}.Files.File{nbfile+1}.TStart=Video.TStart;
s.Animal.Videos.Video{nVideo}.Files.File{nbfile+1}.Duration=round(Video.Duration*1000);
s.Animal.Videos.Video{nVideo}.Files.File{nbfile+1}.SamplingRate=Video.FrameRate;


