function [IdVid,Frame]=GetVideoPosition(CurrExp,CurrBin,handles,params,Timefromstart,CurrVid)

%get start time of all video file and the current abs time
TstartCurrBin=params.FileInfo(CurrExp).BinFiles(CurrBin).TStart*24*3600;
Tcurr=TstartCurrBin+Timefromstart;

allTStartVid=[params.FileInfo(CurrExp).VideosFiles(CurrVid).Files(:).TStart]*24*3600;
allDurationVid=[params.FileInfo(CurrExp).VideosFiles(CurrVid).Files(:).Duration];
if max(allDurationVid(:))==0 %no duration specified inthe exp file
    IdVid=find(Tcurr>=allTStartVid,1,'first');
else
    IdVid=find(Tcurr>=allTStartVid & Tcurr<=allTStartVid+allDurationVid,1,'first');
end
if isempty(IdVid)==0
    %Fps=params.FileInfo(CurrExp).VideosFiles(CurrVid).FrameRate;
    %get the fps from the file
    VideoName=fullfile(params.FileInfo(CurrExp).VideosFiles(CurrVid).Files(IdVid).Dir,params.FileInfo(CurrExp).VideosFiles(CurrVid).Files(IdVid).FileName);
   % filename=params.FileInfo(CurrExp).VideosFiles(CurrVid)
   
   if isempty(params.FileInfo(CurrExp).VideosFiles(CurrVid).Files(IdVid).TimeStamp)==1
        Info=extractAviInfo(VideoName);
        Fps=Info.Fps;
        CurrVideoFiletime=Tcurr-params.FileInfo(CurrExp).VideosFiles(CurrVid).Files(IdVid).TStart*24*3600;
        Frame=round(CurrVideoFiletime*Fps);
   else
       %TimeAbs=Timefromstart/24/3600+params.FileInfo(CurrExp).BinFiles(1).TStart;
       TimerIm=params.FileInfo(CurrExp).VideosFiles(CurrVid).Files(IdVid).TimeStamp;
       Frame=find((TimerIm-(Tcurr/3600/24))>=0,1,'first');
       
   end
    
    %Fps=params.FileInfo(CurrExp).VideosFiles(CurrVid).Files(IdVid).FrameRate;
   
else
    Fps=[];
    CurrVideoFiletime=[];
    Frame=[];
end   