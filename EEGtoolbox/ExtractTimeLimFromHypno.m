function [Header,Data]=ExtractTimeLimFromHypno(Header,MinEpisodDuration,TimeStart,AnalysisDuration,TransitionTime2Remove, ArtefactDetection,StateCode2Process,Event2Remove,SaveAndDisp)
%this function egenerate time limit from an hypnograme

%MinEpisodDuration in sec is the minimale episod duration
%TimeStart time start in min from the beginning of the exp file
%AnalysisDuration total duration to analyse in min
%TransitionTime2Remove in s time to remove when transition in sec, must be a real value
%ArtefactDetection =1 if the artefact lim are used
%StateCode2Process containt the matrix with the state code to process
%Event2Remove containt the event lim to supress from the hypno
%SaveAndDisp is 1 if data save and display 

if isempty(StateCode2Process)
    StateCode2Process=[Header.Info.State(:).Code];
end

if isempty(TransitionTime2Remove)
    TransitionTime2Remove=zeros(1,length(StateCode2Process));
end

if nargin<=7
    Event2Remove=[]; 
    SaveAndDisp=1;
elseif nargin<=8
    SaveAndDisp=1;
end

Data.Settings.TimeLimFromHypno.MinEpisodDuration=MinEpisodDuration;
Data.Settings.TimeLimFromHypno.TransitionTime2Remove=TransitionTime2Remove;
Data.Settings.TimeLimFromHypno.TimeStart=TimeStart;
Data.Settings.TimeLimFromHypno.AnalysisDuration=AnalysisDuration;
Data.Settings.TimeLimFromHypno.ArtefactDetection=ArtefactDetection;


if isempty(Header)==1
     %load expfile
     Info=loadEXP;
     Header.Info=Info;
end
 

%get the absolute time start from the exp file selected
TimeStartExp=Header.Info.BinFiles(1).TStart; %all Tstart are absolute time in day


 %Load Hypno
%compute the whole hypno for all the exp file , the whole hypno begin
%at the first second of the bin file
% LastBinName=[Header.Info.FilesDir Header.Info.BinFiles(end).FileName];
% TotalBinDuration=(Header.Info.BinFiles(end).TStart-Header.Info.BinFiles(1).TStart)*24*3600+GetBinDuration(LastBinName,Header.Info.NbRecChan)/Header.Info.Fs;
% LastHypnoName=[Header.Info.FilesDir Header.Info.HypnoFiles(end).FileName];
% TotalHypnoDuration=(Header.Info.HypnoFiles(end).TStart-Header.Info.HypnoFiles(1).TStart)*24*3600+GetBinDuration(LastHypnoName,1);
% 
% WholeHypno=zeros(ceil(LastHypnoName),1);
% 
% %evaluate the offsetTime
% 
% %append current hypno file to the whole hypno
% for nHypno=1:length(Header.Info.HypnoFiles);
%     CurrHypnoName=[Header.Info.FilesDir Header.Info.HypnoFiles(nHypno).FileName];
%     RelHypnoStart=floor((Header.Info.HypnoFiles(nHypno).TStart-Header.Info.BinFiles(1).TStart)*24*3600);
%     HypFid=fopen(CurrHypnoName);  
%     CurrHypno=fread(HypFid,'uint16');
%     HypnoDuration=length(CurrHypno);
%     fclose(HypFid);
%     WholeHypno(RelHypnoStart+1:RelHypnoStart+HypnoDuration)=CurrHypno;
% end
Params.FileInfo=Header.Info;
Params.FilesDir=Header.Info.FilesDir;
[WholeHypno,TimeScaleAbs,TimeScaleBin,TimeScaleHypno]=ExtractFullHypno(Params,1);




if isinf(AnalysisDuration)==1
       AnalysisDuration=floor((TimeScaleHypno(end))/60-TimeStart);
end

 %read all data
AbsStart=TimeStart*60;
AbsEnd=TimeStart*60+AnalysisDuration*60;


%update the header output data
Header.TStart=AbsStart;
Header.TEnd=AbsEnd;
Header.Fs=Header.Info.Fs;


%current hypno part
%Hypno=WholeHypno(floor(AbsStart)+1:floor(AbsEnd)+1);
%TimeScaleBin=floor(AbsStart)+1:floor(AbsEnd)+1;
%Header.Hypno.Time=linspace(AbsStart,AbsEnd,length(Hypno));
IdHypno=TimeScaleBin>=AbsStart & TimeScaleBin<=AbsEnd;
Header.Hypno.Time=TimeScaleBin(IdHypno);
Header.Hypno.Data=WholeHypno(IdHypno)';
Header.Hypno.RawData=WholeHypno(IdHypno)';
 
 
 
 %artefact detection

 ArtefactMat=ones(length(Header.Hypno.RawData),1);
 if ArtefactDetection==1
     try
        load TempDir
        [FileName,DefaultDir,FilterIndex]=uigetfile({'*.mat'},'Artefact File','multiselect','off',DefaultDir);
    catch
        [FileName,DefaultDir,FilterIndex]=uigetfile({'*.mat'},'Artefact File','multiselect','off','c:\');
    end

    if FilterIndex==1
        load([DefaultDir FileName],'Artef');
    end
    if exist('Artef','var')==0 || FilterIndex~=1 %if artef not load
        %extract artefact
        Artef=ArtefactLim(TimeStart,AnalysisDuration,[],Info);
    end
    
     %add artefact in the hypno
   

    for nartef=1:length(Artef.ArtefTimeDeb)
        IdDebArtef=find(Header.Hypno.Time>=Artef.ArtefTimeDeb(nartef),1,'first');
        IdEndArtef=find(Header.Hypno.Time>Artef.ArtefTimeEnd(nartef),1,'first');

        ArtefactMat(IdDebArtef:IdEndArtef)=0;

    end  
    
    

 end
 
 
 %if some oscillations should be remove
  Hypno2=ones(length(Header.Hypno.RawData),1);
 if isempty(Event2Remove)==0
     
     for nEvent2Remove=1:length(Event2Remove.Deb)
        IdDebNew=find(Header.Hypno.Time>=Event2Remove.Deb(nEvent2Remove),1,'first');
        IdEndNew=find(Header.Hypno.Time>Event2Remove.End(nEvent2Remove),1,'first');

        Hypno2(IdDebNew:IdEndNew)=0;

    end  
     
     
 end
 
 
NewHypno=NaN(length(Header.Hypno.RawData),1);
 Str='';
 for nState=1:length(StateCode2Process)
    CurrStatePos=find([Header.Info.State(:).Code]==StateCode2Process(nState));
    if isempty(CurrStatePos)==0
         CurrHypnoState=Header.Hypno.RawData==Header.Info.State(CurrStatePos).Code;
         Str=[Str Header.Info.State(CurrStatePos).Label '-'];
        %remove transition bouts from the full raw hypno
        hypnoWithZeros=[zeros(TransitionTime2Remove(nState)+2,1);CurrHypnoState.*ArtefactMat.*Hypno2;zeros(TransitionTime2Remove(nState)+1,1)];%add zeros at the begining and end of the file to avoid error when computing the detection
        DiffId=find(diff(hypnoWithZeros)~=0);
        CurrHypno=[];
        if TransitionTime2Remove(nState)>0
            for n=1:length(DiffId)
                hypnoWithZeros(DiffId(n)-TransitionTime2Remove(nState):DiffId(n)+TransitionTime2Remove(nState)-1)=0;
            end
            CurrHypno=hypnoWithZeros(TransitionTime2Remove(nState)+2:end-TransitionTime2Remove(nState)-2);%remove the external zero added for transition detection
            clear hypnoWithZeros;
            %Header.Hypno.Time=linspace(AbsStart,AbsEnd,length(CurrHypno));
        else
            CurrHypno=CurrHypnoState;
        end
        NewHypno(CurrHypno==1)=CurrHypno(CurrHypno==1).*Header.Info.State(CurrStatePos).Code;
        Header.Hypno.Data=NewHypno;
     



         %for each state

         id=(CurrHypno.*ArtefactMat.*Hypno2);
        %  Header.Hypno.Data(logical(id))=NaN;

         %extract artefact start and end times
        StateEpisodStart{nState}=Header.Hypno.Time(diff([0 id'])==1);
        StateEpisodEnd{nState}=Header.Hypno.Time(diff([id' 0])==-1)+1;
        StateEpisodDuration{nState}=StateEpisodEnd{nState}-StateEpisodStart{nState};



         %extract for each state all the episod which last at least MinEpisodDuration


        StateStartEpisodOk{nState}=StateEpisodStart{nState}(StateEpisodDuration{nState}>MinEpisodDuration);
        StateEndEpisodOk{nState}=StateEpisodEnd{nState}(StateEpisodDuration{nState}>MinEpisodDuration);
        StateDurationEpisodOk{nState}=StateEpisodDuration{nState}(StateEpisodDuration{nState}>MinEpisodDuration);
        StateHypno{nState}=ones(1,length(StateStartEpisodOk{nState})).*Header.Info.State(CurrStatePos).Code;
    end
 end
 Data.Event.Deb=[];
 Data.Event.End=[];
 Data.Event.Duration=[];
 Data.Event.Hypno=[];
%  %extract the shorter episod of each state
%  for n=1:4
%      [v,id]=min(StateDurationEpisodOk{n});
%      if isempty(id)==0
%          Data.Event.Deb=cat(1,Data.Event.Deb,StateStartEpisodOk{n}(id));
%          Data.Event.End=cat(1,Data.Event.End,StateEndEpisodOk{n}(id));
% 
%          Data.Event.Duration=cat(1,Data.Event.Duration,StateDurationEpisodOk{n}(id));
%          Data.Hypno=cat(1,Data.Hypno,n);
%      end
% 
%  end
 

 
% concatenate all episod
 [Data.Event.Deb,idorder]=sort([StateStartEpisodOk{:}]);
 Data.Event.End=sort([StateEndEpisodOk{:}]);
 AllDuration=[StateDurationEpisodOk{:}];
 AllState=[StateHypno{:}];
 Data.Event.Duration=AllDuration(idorder);
 Data.Event.Hypno=AllState(idorder);
 Data.Label=Str;
 
 if SaveAndDisp==1
    CurrDir=Header.Info.FilesDir;
    mkdir(CurrDir,'\Computed Data\EventLim\');

    filename=fullfile(CurrDir,'\Computed Data\EventLim\',Header.Info.ExpFileName(1:end-4));

    Time=linspace(0,(length(Header.Hypno.RawData)-1)./60,length(Header.Hypno.RawData))+TimeStart;

    figure;h(1)=subplot(3,1,1);plot(Time,Header.Hypno.RawData);title(h(1),'raw hypno')
    xlabel('Time (Min)')
    h(2)=subplot(3,1,2);
    plot(Time,Header.Hypno.RawData.*ArtefactMat.*Hypno2);title(h(2),'hypno without artefact and without exluded Oscillations');
    h(3)=subplot(3,1,3);plot(Time,Header.Hypno.Data);title(h(3),'Final Hypno')

    linkaxes(h,'x');

    save(sprintf('%s EventLim From Hypno.mat',filename),'Data','Header');
 end
 
 
 
 

 
 
 