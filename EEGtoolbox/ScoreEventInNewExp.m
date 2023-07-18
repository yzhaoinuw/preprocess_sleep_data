function [Info]=ScoreEventInNewExp(OutName, PrefState,NewCode,ExpFileName,DataEvent)
'load the original expfile'

Info=loadEXP(ExpFileName,'No');
Directory=Info.FilesDir;
%duplicate it with the new new'

%Duplicate the hypno and EXP;
OldExpFilenameExt=Info.ExpFileName;
NewEXPFilename=fullfile(Directory,[OutName '.exp']);

[~,NewEXPFilename,~] = fileparts(NewEXPFilename);
NewInfo=Info;
NewEXPFilenameExp=[NewEXPFilename '.exp'];

formatIn = 'yyyy-mm-dd_HH-MM-SS-FFF';


%load old exp struct
[ s ] = Xml2struct(fullfile(Directory,OldExpFilenameExt));

%update the fileorigin in the exp struct
s.Animal.FileOrigin.Text=Directory;

PrefixeFolder=['Hypno_' NewEXPFilename];

mkdir(fullfile(Directory,PrefixeFolder));


for nHyp=1:length(NewInfo.HypnoFiles)
    %update the structure of the nex exp file
    NewInfo.HypnoFiles(nHyp).FileName=fullfile(PrefixeFolder,[NewEXPFilename '_' Info.HypnoFiles(nHyp).FileName(end-length(formatIn)-1:end)]); 
    
    if exist(fullfile(Directory,NewInfo.HypnoFiles(nHyp).FileName))==1
        Ans=questdlg('File already exist, would you like to replace it?');
        if strcmp(Ans,'Yes')    
            %duplicate hypnofile in a new file
            copyfile(fullfile(Info.HypnoFiles(nHyp).Dir,Info.HypnoFiles(nHyp).FileName),...
                fullfile(Directory,NewInfo.HypnoFiles(nHyp).FileName));
        end
    else
        
         copyfile(fullfile(Info.HypnoFiles(nHyp).Dir,Info.HypnoFiles(nHyp).FileName),...
                fullfile(Directory,NewInfo.HypnoFiles(nHyp).FileName));
    end

    %update the new struct
    s.Animal.Hypnogram.Files.File{nHyp}.FileName.Text=NewInfo.HypnoFiles(nHyp).FileName;
   
end


%duplicate the matfile with the new name    
copyfile(fullfile(Directory,[OldExpFilenameExt(1:end-3) 'mat']),...
        fullfile(Directory,[NewEXPFilename '.mat']));
    

%save the new EXP;
XMLCurrFileName=[NewEXPFilename '.xml'];
struct2xml(s,fullfile(Directory,XMLCurrFileName));
copyfile(fullfile(Directory,XMLCurrFileName),fullfile(Directory,NewEXPFilenameExp),'f');
delete(fullfile(Directory,XMLCurrFileName));

%reload the newexpfile;
Info=[];
MaxRangeType='No';
clear Info NewInfo;
Info=loadEXP(fullfile(Directory,NewEXPFilenameExp),MaxRangeType);


%load the event 
if isempty(DataEvent)==1
    [f,MatDir]=uigetfile('* _EventLim.mat','pick the Event file',[Directory]);
    betaLimFile=[MatDir f];
    load(betaLimFile,'Data');
    DataEvent=Data;
end

%score the event

%sort the state

AllEventState=DataEvent.Event.Hypno;
TimeDeb=DataEvent.Event.Deb;
TimeEnd=DataEvent.Event.End;
if isempty(PrefState)==0
    IfOk=ismember(AllEventState,PrefState);
else
    IfOk=true(length(AllEventState),1);
end

for n=1:length(TimeDeb)
    
    if IfOk(n)==1
%         DebHypnoRef=floor(TimeDeb(n));
%         EndHypnoRef=ceil(TimeEnd(n));
        ReplaceStateInHBetweenLim(Info,NewCode,TimeDeb(n),TimeEnd(n));
    end

end

'done'

