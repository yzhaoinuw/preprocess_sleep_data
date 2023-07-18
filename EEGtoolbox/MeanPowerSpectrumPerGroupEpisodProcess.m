function MeanPowerSpectrumPerGroupEpisodProcess(Analyse,Settings,SaveTF)

% Exemple of input parameters
% %PowerSpectrum Setting
% Settings.MTT.Movingwin=[3 0.5];%win and step
% Settings.MTT.BandWidth=2;
% Settings.MTT.Taper=5;
% 
% Settings.MTT.params.trialave = 1;
% Settings.MTT.params.Fs = [];%will be filled after
% Settings.MTT.params.tapers = [Settings.MTT.BandWidth Settings.MTT.Movingwin(1) 2*Settings.MTT.BandWidth*Settings.MTT.Movingwin(1)-Settings.MTT.Taper];
% Settings.MTT.params.fpass = [0.5 45];
% Settings.MTT.params.pad = 1;
% 
% Settings.Filt=[];
% Settings.ApplyICA=0;
% Settings.UseWhitening=0;
% Settings.ARmodel=[];
% Settings.EpochDuration=5;
% Settings.FigSave=0;
% Settings.LightOn=[];
% Settings.ProcessDataPerXMin=[];
% Settings.Artef=[];
% 
% Settings.Band2Analyse={[1 4] [4 9]};
% 
% Analyse(1).Mice={'S1' 'S2' 'S3' 'S4' 'S6' 'S7' 'S8' 'S9' 'S10' 'S11' 'S12'};
% Analyse(1).Baseline='BL (10H-16H)-SWS';
% Analyse(1).Manip='APSD (10H-16H)-SWS';
% Analyse(1).Duration=360;
% Analyse(1).TransitionTime2Remove=10;
% 
% Analyse(2).Mice={'S1' 'S2' 'S3' 'S4' 'S6' 'S7' 'S8' 'S9' 'S10' 'S11' 'S12'};
% Analyse(2).Baseline='BL (16H-20H)-SWS';
% Analyse(2).Manip='PSR (16H-20H)-SWS';
% Analyse(2).Duration=240;
% Analyse(2).TransitionTime2Remove=10;
% 
% Analyse(3).Mice={'S1' 'S2' 'S3' 'S4' 'S6' 'S7' 'S8' 'S9' 'S10' 'S11' 'S12'};
% Analyse(3).Baseline='BL (16H-20H)-PS';
% Analyse(3).Manip='PSR (16H-20H)-PS';
% Analyse(3).Duration=240;
% Analyse(3).TransitionTime2Remove=10;
% 
% %Settings extraction states
% Settings.MinEpisodDuration=0;%0s  %5s in min = 0.083333333333333
% Settings.ArtefactDetection=0;
% 
% SaveTF=0


close all;



%pick the xls file that cntain the infomration to process
try     load TempDir 
    [f,DefaultDir]=uigetfile({'*.xls;*.xlsx','xls files'},'pick the xls file',DefaultDir);
catch
    [f,DefaultDir]=uigetfile({'*.xls;*.xlsx','xls files'},'pick the xls file');
end
save TempDir DefaultDir;


XlsFileName=fullfile(DefaultDir,f);

[~,~,MAT]=xlsread(XlsFileName,'Feuil1');

AllGroupList=MAT(2:end,3);
AllAnimalList=MAT(2:end,2);
AllSubDir=MAT(2:end,4);
AllDExpFileName=MAT(2:end,5);
AllT0=[MAT{2:end,6}];
%AllTEnd=[MAT{2:end,7}]; not use
AllChannel=[MAT{2:end,8}];
AllState2Analyse=cat(1,MAT{2:end,10});

Maniptype=unique(AllGroupList);

AnimalList=unique(AllAnimalList);



%create empty XlsREsultFile
NameXlsOut=['Result' f];


PwBL={};
PwManip={};
PwNorm={};
MeanBandPwBL={};
MeanBandPwManip={};
F={};
PWNormTab=[];
DataAllLim={};


for nAnalyse=1:length(Analyse)  %for each analyse
    
    
    
    PwBL={};
    PwManip={};
    PwNorm={};
    MeanBandPwBL={['Mean Band PW/ ' Analyse(nAnalyse).Baseline]};
    MeanBandPwManip={['Mean Band PW/ ' Analyse(nAnalyse).Manip]};
    
    for nband=1:length(Settings.Band2Analyse)
        
        MeanBandPwBL{1,nband+1}=num2str(Settings.Band2Analyse{nband});
        MeanBandPwManip{1,nband+1}=num2str(Settings.Band2Analyse{nband});
        
    end
    
    F={};
    PWNormTab=[];
    
    
    
   disp(sprintf('///////// Manip: %s',Analyse(nAnalyse).Manip))

   
   %initfigure
   Fig(nAnalyse)=figure('name',Analyse(nAnalyse).Manip,'units','pixels','position',[100 100 1000 1000]);
   
   clear Ha;
   
   Ha(1,1)=subplot(2,1,1);
   Ha(1,2)=subplot(2,1,2);
   linkaxes(Ha,'x');
   
   Mousecolor=hsv(length(Analyse(nAnalyse).Mice));
   MouseOk=[];
   
   
   clear Duration TransitionTime2Remove;
   Duration=Analyse(nAnalyse).Duration;
   TransitionTime2Remove=Analyse(nAnalyse).TransitionTime2Remove;
   
   
    %extract timefrequency for baseline
    for nmouse=1:length(Analyse(nAnalyse).Mice)
        disp(sprintf('Mouse: %s',Analyse(nAnalyse).Mice{nmouse}))
        idmousebaseline=find(ismember(AllAnimalList,Analyse(nAnalyse).Mice{nmouse}) & ismember(AllGroupList,Analyse(nAnalyse).Baseline));
        idmousemanip=find(ismember(AllAnimalList,Analyse(nAnalyse).Mice{nmouse}) & ismember(AllGroupList,Analyse(nAnalyse).Manip));
        MouseOk=[MouseOk idmousemanip];
        if length(idmousebaseline)==1 && length(idmousemanip)==1
            %load file
            currfilename=fullfile(DefaultDir,AllSubDir{idmousebaseline},AllDExpFileName{idmousebaseline});
            HeaderBaseline.Info=loadEXP(currfilename,'No');
            
            
           
            StateCode2Process=AllState2Analyse(idmousebaseline);
            if isnumeric(StateCode2Process)==0
                StateCode2Process=str2double(StateCode2Process);
            end
            
            CodeLabelList={HeaderBaseline.Info.State(:).Label};
            CodeColorList=cat(1,HeaderBaseline.Info.State(:).Color);
            StateName=CodeLabelList{ismember([HeaderBaseline.Info.State.Code],StateCode2Process)};
            StateColor=CodeColorList(ismember([HeaderBaseline.Info.State.Code],StateCode2Process),:);
            %
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %extact Hypno limit from baseline
            T0=AllT0(idmousebaseline);
            %Duration=AllTEnd(idmousebaseline)-T0;
            
            clear TF F DataLimState;
           
            [~,DataLimState]=ExtractTimeLimFromHypno(HeaderBaseline,Settings.MinEpisodDuration,T0,Duration,TransitionTime2Remove, Settings.ArtefactDetection,StateCode2Process,[],0);
            DataAllLim{nAnalyse,nmouse}.Baseline=DataLimState;
            
            ChanNum=AllChannel(idmousemanip);
            [TF,f]=EpisodTimeFrequency(DataLimState.Event,HeaderBaseline.Info,ChanNum,Settings.MTT);
 
            Analyse(nAnalyse).F{nmouse}=f;
            
            TFBaseline{nmouse}.S=TF;
            
            if SaveTF==1
                Analyse(nAnalyse).TFBaseline{nmouse}=TFBaseline{nmouse};
                Analyse(nAnalyse).HypnoLimBaseline{nmouse}=DataLimState;
            end
            
            F=[{'F (Hz)'} num2cell(Analyse(nAnalyse).F{nmouse})];

            Analyse(nAnalyse).MeanPWBL{nmouse}=mean(TFBaseline{nmouse}.S,1);
            NbLine=size(PwBL,1);
            PwBL{NbLine+1,1}={['PW-' Analyse(nAnalyse).Mice{nmouse}]};
            PW=num2cell(Analyse(nAnalyse).MeanPWBL{nmouse});
            PwBL(NbLine+1,2:length(PW)+1)=PW;

            MeanBandPwBL{NbLine+2,1}=[Analyse(nAnalyse).Mice{nmouse}];
            
            for nband=1:length(Settings.Band2Analyse)
                idband=f>=Settings.Band2Analyse{nband}(1) & f<=Settings.Band2Analyse{nband}(2) ;
                Analyse(nAnalyse).MeanBandPowerBL(nmouse,nband)=mean(Analyse(nAnalyse).MeanPWBL{nmouse}(idband));
                MeanBandPwBL{NbLine+2,nband+1}=Analyse(nAnalyse).MeanBandPowerBL(nmouse,nband);
            end 

			line(Ha(1,1),'xdata',Analyse(nAnalyse).F{nmouse},'ydata',10*log10(Analyse(nAnalyse).MeanPWBL{nmouse}),'color',Mousecolor(nmouse,:),'linestyle','--');
           % line(Ha(1,1),'xdata',Analyse(nAnalyse).F{nmouse},'ydata',Analyse(nAnalyse).MeanPWBL{nmouse},'color',Mousecolor(nmouse,:),'linestyle','--');
            grid(Ha(1,1),'on')

           

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %manip
            %load file
            currfilenameManip=fullfile(DefaultDir,AllSubDir{idmousemanip},AllDExpFileName{idmousemanip});
            HeaderManip.Info=loadEXP(currfilenameManip,'No');
            
            StateCode2Process=AllState2Analyse(idmousemanip);
            if isnumeric(StateCode2Process)==0
                StateCode2Process=str2double(StateCode2Process);
            end
            
            CodeLabelList={HeaderManip.Info.State(:).Label};
            CodeColorList=cat(1,HeaderManip.Info.State(:).Color);
            StateName=CodeLabelList{ismember([HeaderManip.Info.State.Code],StateCode2Process)};
            StateColor=CodeColorList(ismember([HeaderManip.Info.State.Code],StateCode2Process),:);
            
            %extact Hypno limit from manip
            
            T0=AllT0(idmousemanip);
            %Duration=AllTEnd(idmousebaseline)-T0;
            %Duration=6*60;
            clear TF DataLimState;
            [~,DataLimState]=ExtractTimeLimFromHypno(HeaderManip,Settings.MinEpisodDuration,T0,Duration,TransitionTime2Remove, Settings.ArtefactDetection,StateCode2Process,[],0);
            DataAllLim{nAnalyse,nmouse}.Manip=DataLimState;
     
            [TF,f]=EpisodTimeFrequency(DataLimState.Event,HeaderManip.Info,ChanNum,Settings.MTT);

            TFManip{nmouse}.S=TF;
            
            
            
            if SaveTF==1
                Analyse(nAnalyse).TFManip{nmouse}=TFManip{nmouse};
                Analyse(nAnalyse).HypnoLimManip{nmouse}=DataLimState;
            end

            Analyse(nAnalyse).MeanPWManip{nmouse}=mean(TF,1);

            NbLine=size(PwManip,1);
            PwManip(NbLine+1,1)={['PW-' Analyse(nAnalyse).Mice{nmouse}]};
            PW=num2cell(Analyse(nAnalyse).MeanPWManip{nmouse});

            MeanBandPwManip{NbLine+2,1}=[Analyse(nAnalyse).Mice{nmouse}];

             for nband=1:length(Settings.Band2Analyse)
                idband=f>=Settings.Band2Analyse{nband}(1) & f<=Settings.Band2Analyse{nband}(2);
                Analyse(nAnalyse).MeanBandPowerManip(nmouse,nband)=mean(Analyse(nAnalyse).MeanPWManip{nmouse}(idband));
                MeanBandPwManip{NbLine+2,nband+1}=Analyse(nAnalyse).MeanBandPowerManip(nmouse,nband);
             end
             
              L(nmouse)=line(Ha(1),'xdata',Analyse(nAnalyse).F{nmouse},'ydata',10*log10(Analyse(nAnalyse).MeanPWManip{nmouse}),'color',Mousecolor(nmouse,:),'linestyle','-');
			  %L(nmouse)=line(Ha(1),'xdata',Analyse(nAnalyse).F{nmouse},'ydata',Analyse(nAnalyse).MeanPWManip{nmouse},'color',Mousecolor(nmouse,:),'linestyle','-');
            
              %normalisation from baseline
             Analyse(nAnalyse).MeanNormPWManip{nmouse}=(Analyse(nAnalyse).MeanPWManip{nmouse}-Analyse(nAnalyse).MeanPWBL{nmouse})./Analyse(nAnalyse).MeanPWBL{nmouse}.*100;
             PwNorm(NbLine+1,1)={['nPW-' Analyse(nAnalyse).Mice{nmouse}]};
             PwNorm(NbLine+1,2:length(PW)+1)=num2cell(Analyse(nAnalyse).MeanNormPWManip{nmouse});
             PWNormTab(NbLine+1,1:length(PW))=Analyse(nAnalyse).MeanNormPWManip{nmouse};


             set(Ha(1),'xlim',[0 Analyse(nAnalyse).F{nmouse}(end)]);
             set(Ha(2),'xlim',[0 Analyse(nAnalyse).F{nmouse}(end)]);
             grid(Ha(2),'on')
      
     
        else
            disp(sprintf('mouse %s,data missing',Analyse(nAnalyse).Mice{nmouse}))
        end
    
        drawnow;
    
    end
   
    PwNorm =cat(1,F,PwNorm);
        
    MeanNPW=mean(PWNormTab(:,:),1);
    StdNPW=std(PWNormTab(:,:),1,1);
    SemNPW=StdNPW./sqrt(size(PWNormTab,1));
    
    StateColor=[0 0 0];
    
    line(Ha(2),'xdata',Analyse(nAnalyse).F{1},'ydata',MeanNPW,'color',StateColor,'linestyle','-');

    patch(Ha(2),'XData',[Analyse(nAnalyse).F{1} fliplr(Analyse(nAnalyse).F{1})],'YData',[MeanNPW+SemNPW fliplr(MeanNPW-SemNPW)],'Facecolor',StateColor,'FaceAlpha',0.1,'LineStyle','none')
    set(Ha(2),'ylim',[-200 200]);

   % errorbar(Analyse(nAnalyse).F{1},MeanNPW,SemNPW,'color',StateColor(nstate,:));
    box(Ha(2),'off')
    box(Ha(1),'off')

    xlabel(Ha(1),'Frequency (Hz)');
    ylabel(Ha(1),sprintf(' PWS (dB)'));
    title(Ha(1),sprintf('%s RAW PWS %s',Analyse(nAnalyse).Manip,StateName));
    legend(L,AllAnimalList(MouseOk));

    xlabel(Ha(2),'Frequency (Hz)');
    ylabel(Ha(2),sprintf('%% of change from baseline'));
    title(Ha(2),sprintf('%s Norm PWS %s',Analyse(nAnalyse).Manip,StateName));
        

    [Directory,FileName]=fileparts(XlsFileName);
    FigFileName=fullfile(Directory,[FileName '-' Analyse(nAnalyse).Manip '.fig']) ;
    BmpFileName=fullfile(Directory,[FileName '-' Analyse(nAnalyse).Manip '.bmp']); 
    pdfFileName=fullfile(Directory,[FileName '-' Analyse(nAnalyse).Manip '.pdf']); 
    MatFileName=fullfile(Directory,[FileName '-' Analyse(nAnalyse).Manip '.mat']); 
    
    
    hgsave(Fig(nAnalyse), FigFileName, '-v7.3');   
    set(Fig(nAnalyse),'PaperPositionMode','Auto');
    print(Fig(nAnalyse),BmpFileName,'-dbmp','-r100');
    print(Fig(nAnalyse),pdfFileName,'-dpdf','-r200','-bestfit');
    save(MatFileName,'Analyse','DataAllLim');
    
    xlswrite(XlsFileName,PwNorm,[ Analyse(nAnalyse).Manip '- Norm PW']);
    xlswrite(XlsFileName,[MeanBandPwBL MeanBandPwManip],[ Analyse(nAnalyse).Manip '- Band PW']);
    
    

    
end


