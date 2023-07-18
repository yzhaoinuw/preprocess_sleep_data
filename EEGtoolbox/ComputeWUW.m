function [MatFileName,Header,HB,Raw,TF,Artef,Filt,BandPW,PeakFqcy,RatioBandPW,Peak,XCor,HilbertReal,Cohe,ARModel]=ComputeWUW(OutNamePref,T0,TotDuration,ProcessInfo,Pass,Info,Artef,FigSave,BoutSize,NBTranstionWin2Remove,ProcessDataPerXMin,Filt,MTT,LightOn,MatFileNamePass1,ARModel,MainOutDir,PerSat2Remove)
%warning off all;
%OutNamePref is a prefix for the output files
%To is the time to start the analyse in min from tstart form the exp selected
%TotDuration is the duration to analyse in min
%Pass is 1 if first analyse (time frequency)
%Pass is 2 is second analyse compute  band sum, ratio band on the TF or peak detection one the cross corrn TF or COV should exist
%INfo contain the expfile (obtain from the function loadEXP, could be the value should be 1 if not 0
%BoutSize=5; %epoch duration  in s
%NBTranstionWin2Remove=0; %number of epoch to remove   
%ProcessDataPerXMin=5; %process data per step of ProcessDataPerXMin if empty do not cut the nalaysis
%Filt is ta structure that containt the different filter to use, could be empty if no filter is used
%MTT is ta structure that containt the different Multitaper settings to use, could be empty if no filter is used
%LightOn could be empty if not should be a mtrix with the start hour and end hour of light on
%MatFileNamePass1 the full filename of the pass1 to load
%ARModel is the model to use for whiteneing could be empty
%CurrDir=[cd '\'];
    %warning off all;

   % fclose all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load expfile
    if exist('Info')==0 || (exist('Info')==1  && isempty(Info)==1);
        Info=loadEXP;
    end
  CurrDir=Info.FilesDir;  
  
  %create folder %s\\Computed Data
  mkdir(sprintf('%s\\Computed Data',CurrDir));
  
  if nargin<=13
      LightOn=[];
      MatFileNamePass1='';
      ARModel=[];
      PerSat2Remove=[];
  elseif nargin==14
      MatFileNamePass1='';
      ARModel=[];
       PerSat2Remove=[];
  elseif nargin==15
      ARModel=[];
      MainOutDir=[];
       PerSat2Remove=[];
      
  elseif nargin==16
      MainOutDir=[];
       PerSat2Remove=[];
  elseif nargin==17
       PerSat2Remove=[];
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %define generic parameter
    
    %export parameters
    %ExportXLS=1; %ExportXLS is 1 if data are exported in xls 0 if not 2 if exported in xlsx
    MatSave=1; %if data must be saved in a mat file

    %process is define a describe here
    %first cell is the channel if 2 value are in a line vector (ex: '[2 3]') means
    %a difference if 2 value are in a row vector (ex: '[2;3]') means 2 chan process for covariance for ex
    
    %second cell is the kind of parameter to evaluate
    %   1 Compute TimeFrequency with MTT
    %   2 is heart rate  
    %   3 raw signal
    %   4 Peak frequency in specfic band
    %   5 mean Band power
    %   6 pattern correlation
    %   7 Band ratio
    %   8 Peak detection
    %   9 Peak detection from covariance between 2 chan
    %   10 Hilbert real part
    %   11 Band ratio band sup band inf
    %   12 coherence (phase ampl cross)


    %the third cell is the option for each parameter
    %   1: MTT parameters number
    %   2: [minpeakdist (in s) minpeakvalue]
    %   3: [] 
    %   4: [fmin1 fmax1]; use only if Time frequency have been calculated
    %   5: [fmin1 fmax1]; use only if Time frequency have been calculated
    %   6: mat file which contain the reference pattern in a matrix named
    %   Pattern with a [m,1] size the mat file should be in the same  directory of the exp file pattern is in uint16 value
    %   7: [fmin1 fmax1;fmin2 fmax2] use only if Time frequency have been calculated
    %   8: [MinPeakDist MinPeakThreshold] MinPeakDist is the minimum time in s between two peaks,MinPeakThreshold if thresh is negative find negative peak, MinPeakThreshold is in number of std from the mean
    %   9: [Win Step MinPeakDist MinPeakThreshold]
            %MinPeakDist is the minimum time in s between two peaks,
            %MinPeakThreshold if thresh is negative find negative peak MinPeakThreshold is in number of std from the mean
            %Win is the win size for covariance evaluation in s
            %Step is the step for the process in s
            %id pass =1 MinPeakDist and MinPeakThreshold are unused, set these parameters to NaN
            %if pass 2 just peak is extracted from covariance already done just set MinPeakDist and MinPeakThreshold
    %   10:[fmin1 fmax1] band filtering
    %   11:[fmin1 fmax1 fmin2 fmax2 fmin3 fmax3] f1 /mean(f2;f3)use only if Time frequency have been calculated
    %   12 :MTT parameters number

    %the fourth cell is the decimate factor %if empty mean no subsampling

    %the fifth cell is the number of filter to use, if empty no filter else
    %should be 1 or 2 depending on filter designed above or [] if no filter

    %the sixth cell is the ylim for display or clim depending on the kind
    %of analysis, if empty limit is set automatically
    
    %the seventh cell is if yes or no artefact removing is used  (base on ICA detection);
    
    %the eighth cell is 1 if yes or no Ica Correction is applyed; this value could be empty for Band power calculation

    %the nine cell is 0 if the signal should not be whitened 1 or not 0
    
    %the tenth cell is empty or 0 if no normalisation on the whole TF 1 if
    %normalisation is on. use only for process 4 5 7 and 11 pass2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% process
    
    %init the number of total output variable
    nbVar=0;
    %init the number of total output varname
    VarName={};


    %evaluate the number of variables and generate variable names
    for nProc=1:size(ProcessInfo,1)
        
        if length(ProcessInfo(nProc,:))==9 || (length(ProcessInfo(nProc,:))==10 && isempty(ProcessInfo(nProc,:))==1)
            ProcessInfo{nProc,10}=0;
        end
        
        
        if length(ProcessInfo{nProc,1})==1 %if no substraction 
            ChName=Info.ChLabel{ProcessInfo{nProc,1}(1)};
        elseif size(ProcessInfo{nProc,1},2)==2 %if substraction 
            ChName=sprintf('%s-%s',Info.ChLabel{ProcessInfo{nProc,1}(1)},Info.ChLabel{ProcessInfo{nProc,1}(2)});
        elseif size(ProcessInfo{nProc,1},1)==2 %if 2 chan for covariance peak detection 
            ChName=sprintf('%s/%s',Info.ChLabel{ProcessInfo{nProc,1}(1)},Info.ChLabel{ProcessInfo{nProc,1}(2)});
        end

        if ProcessInfo{nProc,2}==1 %if band analysis           
            nbVar=nbVar+1;
            VarName{nbVar}=sprintf('%s TF',ChName );
        end

        if ProcessInfo{nProc,2}==2 %if HeartRate
                nbVar=nbVar+1;
                VarName{nbVar}=sprintf('%s HR(BPM)',ChName);
        end
        
        if ProcessInfo{nProc,2}==3 %if raw data
                nbVar=nbVar+1;
                VarName{nbVar}=sprintf('%s',ChName);
        end
        
        if ProcessInfo{nProc,2}==4 %if Frequency peak data
                nbVar=nbVar+1;
                VarName{nbVar}=sprintf('%s Pk %5.1f-%5.1f',ChName,ProcessInfo{nProc,3}(1),ProcessInfo{nProc,3}(2));
        end
        
        if ProcessInfo{nProc,2}==5 %if band Sum
                nbVar=nbVar+1;
                VarName{nbVar}=sprintf('%s PW %4.1f-%4.1f',ChName,ProcessInfo{nProc,3}(1),ProcessInfo{nProc,3}(2));
        end
        
        if ProcessInfo{nProc,2}==6 %if pattern correlation
                nbVar=nbVar+1;
                VarName{nbVar}=sprintf('%s PC',ChName);
        end
        
        if ProcessInfo{nProc,2}==7 %if band power ratio
                nbVar=nbVar+1;
                VarName{nbVar}=sprintf('%s PW %4.1f-%4.1f/%4.1f-%4.1f',ChName,ProcessInfo{nProc,3}(1),ProcessInfo{nProc,3}(2),ProcessInfo{nProc,3}(3),ProcessInfo{nProc,3}(4));
        end
        
        if ProcessInfo{nProc,2}==8 %if Peak detection
                nbVar=nbVar+1;
                VarName{nbVar}=sprintf('%s min Peak Thresh %4.1f',ChName,ProcessInfo{nProc,3}(2));
        end
        
        if ProcessInfo{nProc,2}==9 %if Peak detection from XCorr between 2 chan
                nbVar=nbVar+1;
                VarName{nbVar}=sprintf('%s XCorr Peak  (Win %4.1f, step %4.1f, std x %4.1f)',ChName,ProcessInfo{nProc,3}(1),ProcessInfo{nProc,3}(2),ProcessInfo{nProc,3}(4));
        end
        
        if ProcessInfo{nProc,2}==10 %Hilbert real part
                nbVar=nbVar+1;
                VarName{nbVar}=sprintf('%s Hilb Real %4.1f-%4.1f',ChName,ProcessInfo{nProc,3}(1),ProcessInfo{nProc,3}(2));
        end
        
        if ProcessInfo{nProc,2}==11 % band power ratio with 3 band lim corrected by the frequency
                nbVar=nbVar+1;
                VarName{nbVar}=sprintf('%s (%4.1f-%4.1f)/((%4.1f-%4.1f);(%4.1f-%4.1f))',ChName,ProcessInfo{nProc,3}(1),ProcessInfo{nProc,3}(2),ProcessInfo{nProc,3}(3),ProcessInfo{nProc,3}(4),ProcessInfo{nProc,3}(5),ProcessInfo{nProc,3}(6));
        end
        
        if ProcessInfo{nProc,2}==12 % Coherence
                nbVar=nbVar+1;
                VarName{nbVar}=sprintf('%s Coherence',ChName);
        end
        



    end
    
   


    %get the absolute time start from the exp file selected
    TimeStartExp=Info.BinFiles(1).TStart; %all Tstart are absolute time in day
    
    
    %compute the whole hypno for all the exp file , the whole hypno begin
    %at the first second of the bin file
    
    
    
    LastBinName=fullfile(Info.BinFiles(end).Dir,Info.BinFiles(end).FileName);
    %TotalBinDuration=(Info.BinFiles(end).TStart-Info.BinFiles(1).TStart)*24*3600+GetBinDuration(LastBinName,Info.NbRecChan)/Info.Fs;
    TotalBinDuration=etime(datevec(Info.BinFiles(end).TStart),datevec(Info.BinFiles(1).TStart))+GetBinDuration(LastBinName,Info.NbRecChan)/Info.Fs;

    
   %get fullhypno=
    params.FileInfo=Info;
    [WholeHypno,~,TimeScaleBin,~]=ExtractFullHypno(params,1);
    WholeHypno=WholeHypno';
    TimeScaleBin=TimeScaleBin';
%     
%     
%     
%     WholeHypno=zeros(ceil(TotalBinDuration),1);
% 
%     
%     HypOffset=etime(datevec(Info.HypnoFiles(1).TStart),datevec(Info.BinFiles(1).TStart));
%     %append current hypno file to the whole hypno
%     for nHypno=1:length(Info.HypnoFiles);
%         CurrHypnoName=fullfile(Info.HypnoFiles(nHypno).Dir,Info.HypnoFiles(nHypno).FileName);
%         %RelHypnoStart=floor(etime(datevec(Info.HypnoFiles(nHypno).TStart),datevec(Info.HypnoFiles(1).TStart)));
%         RelHypnoStart=floor(etime(datevec(Info.HypnoFiles(nHypno).TStart),datevec(Info.BinFiles(1).TStart)));
%         HypFid=fopen(CurrHypnoName); 
%         CurrHypno=fread(HypFid,'uint16');
%         HypnoDuration=length(CurrHypno);
%         fclose(HypFid);
%         WholeHypno(RelHypnoStart+1:RelHypnoStart+HypnoDuration)=CurrHypno;
%     end

    %remove transition bouts from the full raw hypno
    hypnoWithZeros=[zeros((NBTranstionWin2Remove+1)*BoutSize+1,1);WholeHypno;zeros((NBTranstionWin2Remove+1)*BoutSize,1)];%add zeros at the begining and end of the file to avoid error when computing the detection
    DiffId=find(diff(hypnoWithZeros)~=0);

    if NBTranstionWin2Remove>0
        for n=1:length(DiffId)
            hypnoWithZeros(DiffId(n)-(NBTranstionWin2Remove)*BoutSize:DiffId(n)+NBTranstionWin2Remove*BoutSize-1)=0;
        end
        NewHypno=hypnoWithZeros((NBTranstionWin2Remove+1)*BoutSize+1:end-(NBTranstionWin2Remove+1)*BoutSize-1);%remove the external zero added for transition detection
    else
        NewHypno=WholeHypno;
    end 
    clear WholeHypno;

   if isinf(TotDuration)==1
       TotDuration=floor((length(NewHypno-1)-T0*60)/60)-10/60;
   end
        
     %if the analyse must be cut
    if isempty(ProcessDataPerXMin)==0
        T0Cut=T0:ProcessDataPerXMin:T0+TotDuration;
        if T0Cut(end)~=T0+TotDuration
            T0Cut=[T0Cut T0+TotDuration]; 
        end
        TotDurationCut=diff(T0Cut);
        T0Cut(end)=[];
    else
        T0Cut=T0;
        TotDurationCut=TotDuration;
    end
    
    


     barre=waitbar(0,sprintf('Keep Cool, %d%% processed',0));
     
     Header=[];
     Header.Process=ProcessInfo;
     Header.Info=Info;
     
    
    %init the out put variable
    TF=[];
    HB=[];
    %Artef=[];
    BandPW=[];
    PeakFqcy=[];
    RatioBandPW=[];
    Peak=[];
    XCor=[];
    HilbertReal=[];
    Raw=[];
    Cohe=[];
     
     %start Analysis per period
    for nCutPeriod=1:length(T0Cut)
        
        
        if Pass==2 %if second pass
            %try to load time frequency
            BandPW=[];
            PeakFqcy=[];
            RatioBandPW=[];
            
            if isempty(MatFileNamePass1)==0 && exist(MatFileNamePass1)==2
                try
                    load(MatFileNamePass1,'TF');
                catch
                    TF=[];
                end
                try
                    load(MatFileNamePass1,'XCor');
                end
                 try
                    load(MatFileNamePass1,'Cohe');
                end
                'mat file loaded!'; 
            else
                
                MatFileName=sprintf('%s\\Computed Data\\%s %s process %06.1f min-%06.1f min Pass1.mat',CurrDir,OutNamePref,Info.ExpFileName(1:end-4),T0Cut(nCutPeriod),T0Cut(nCutPeriod)+TotDurationCut(nCutPeriod));


                if exist(MatFileName);
                    try
                        load(MatFileName,'TF');
                    catch
                        TF=[];
                    end
                    try
                        load(MatFileName,'XCor');
                    end
                     try
                        load(MatFileName,'Cohe');
                    end
                    'mat file loaded!'; 

                else
                    'mat file not loaded pick a file!'
                    if exist([CurrDir '\' OutNamePref])
                        [f,MatDir]=uigetfile('* Pass1.mat','pick the Pass 1',[CurrDir '\' OutNamePref]);
                    else
                        [f,MatDir]=uigetfile('* Pass1.mat','pick the Pass 1',[CurrDir]);
                    end
                    MatPass1file=[MatDir f];
                    try
                        load(MatPass1file,'TF');
                    catch
                        TF=[];
                    end
                    try
                    load(MatPass1file,'XCor');
                    end
                    try
                    load(MatPass1file,'Cohe');
                    end


                end
            end
           

        end


        %read all data
        AbsStart=T0Cut(nCutPeriod)*60;
        AbsEnd=(T0Cut(nCutPeriod)+TotDurationCut(nCutPeriod))*60;
       
        

        %update the header output data
        Header.TStart=AbsStart;
        Header.TEnd=AbsEnd;
        Header.Fs=Info.Fs;
        %current hypno part
        
        IdHypno=TimeScaleBin>=AbsStart & TimeScaleBin<AbsEnd;
        %CurrHypno=NewHypno(IdHypno);
        Header.Hypno.Time=(TimeScaleBin(IdHypno));
        Header.Hypno.Data=NewHypno(IdHypno);
        
        
        %CurrHypno=NewHypno(floor(AbsStart)+1:floor(AbsEnd)+1);

        
%         %here we need to add the offset of the hypno due to the round of
%         %5sec
%         Correction=HypOffset-floor(HypOffset);
%         Header.Hypno.Time=(linspace(AbsStart,AbsEnd,length(CurrHypno))+Correction)';
%         Header.Hypno.Data=CurrHypno;



 



        %artefact evaluation from ICA analysis
        if  max([ProcessInfo{:,7}])==1 %if at least one process use artefact
            
            if exist('Artef','var')==0
                Artef=[];
            end
            if  isempty(Artef)==1
                try
                    load TempDir
                    [FileName,DefaultDir,FilterIndex]=uigetfile({'*.mat'},'Artefact File','multiselect','off',DefaultDir);
                catch
                    [FileName,DefaultDir,FilterIndex]=uigetfile({'*.mat'},'Artefact File','multiselect','off','c:\');
                end

                if FilterIndex==1
                    load([DefaultDir FileName],'Artef');
                end

                if  isempty(Artef)==1 || FilterIndex~=1 %if artef not load
                    %extract artefact
                    Artef=ArtefactLim(T0,TotDuration,[],Info);
                end
                
                
               
                
                
                
                
            end
            
            
         %score artefact (5) the hypno if artefact is apply

         if isempty(Artef.ArtefTimeDeb)==0 && isempty(Artef.ArtefTimeEnd)==0
             for nArtef=1:length(Artef.ArtefTimeDeb)
                 id=find(Header.Hypno.Time>=Artef.ArtefTimeDeb(nArtef) & Header.Hypno.Time<=Artef.ArtefTimeEnd(nArtef));
                 Header.Hypno.Data(id)=5;
             end
         end
            
          

            
            
        else
            
            %init artefact time
            Artef.ArtefTimeDeb=[];   
            Artef.ArtefTimeEnd=[]; 

        end 
        
        
        
       



        %init the figure which containt the processes
        if FigSave==1
            hf(nCutPeriod)=figure;
            set(hf(nCutPeriod),'units','normalized');
            set(hf(nCutPeriod),'units','normalized','position',[0 0 1 1]);

            Xo=0.05;
            Yo=0.02;
            W=0.93;
            Sp=0.005;
            H=(1-(2*Yo))./(size(ProcessInfo,1)+1)-Sp;
            h(1,nCutPeriod)=subplot('position',[Xo,1-Yo-H-Sp,W,H]);
            plot(Header.Hypno.Time./60,Header.Hypno.Data);ylabel ('Hypno');
            grid on;
            set(h(1,nCutPeriod),'Ylim',[min([Header.Info.State.Code])-1 max([Header.Info.State.Code])+1],'xticklabel',[]);

        end

        %start each process
        for nprocess=1:size(ProcessInfo,1)
 
            
            try
                
                %exctract the uint16 raw data from bin file
                numchan=ProcessInfo{nprocess,1};
                
                if size(numchan,2)==1% if no substraction
                    if size(numchan,1)==1
                        numchan(2)=numchan(1);
                    else %if convolution without substraction
                        numchan(:,2)=numchan(:,1);
                    end
                end
                
                
                 if Pass==2 && isfield(TF.Chan(numchan(1,1),numchan(1,2)),'S')==1 && isempty(TF.Chan(numchan(1,1),numchan(1,2)).S)==0
                     if ismember(ProcessInfo{nProc,2},[4 5 7 11])==1 && ProcessInfo{nProc,10}==1 %if normalisation ==1
                    
                        MeanPS=nanmean(TF.Chan(numchan(1,1),numchan(1,2)).S,1);
                        TF.Chan(numchan(1,1),numchan(1,2)).S=TF.Chan(numchan(1,1),numchan(1,2)).S./ ...
                            repmat(MeanPS,size(TF.Chan(numchan(1,1),numchan(1,2)).S,1),1);
                        
                     end

                end
     
                
                
                if  ((ProcessInfo{nprocess,2}~=4 && ProcessInfo{nprocess,2}~=5 && ProcessInfo{nprocess,2}~=7 && ProcessInfo{nprocess,2}~=9 && ProcessInfo{nprocess,2}~=11) && Pass==2) ||  Pass==1;%extract data from bin file if not working on TF or pass 2

              

                    if numchan(1,1)==numchan(1,2) %if no substraction
                       

                        if ProcessInfo{nprocess,8}==1 %if ICA correction
                            AllBinData=ExtractContinuousData(Info.FilesDir,Info,[],AbsStart,AbsEnd,PerSat2Remove);
                            AllBinDataICA=AllBinData;
                            
                            Mat_ICA_corr=Info.Artef.ICAMat(:,:,2);

                            % %get the chan of interest and chan 2 correct
                            ChanOfInterest=eval(Info.Artef.ChanOfInterest); 
                            ChanToCorrect=eval(Info.Artef.ChanToCorrect);

                            AllBinDataICA(ChanOfInterest,:) = (pinv(Mat_ICA_corr)*(Mat_ICA_corr* AllBinData(ChanOfInterest,:)));

                            %get the value not corrected
                            AllBinData(ChanToCorrect,:)=AllBinDataICA(ChanToCorrect,:);

                            BinData=AllBinData(numchan(1,1),:);
                            
                            clear AllBinDataICA AllBinData;
                        else
                            BinData=ExtractContinuousData(Info.FilesDir,Info,numchan(1,1),AbsStart,AbsEnd,PerSat2Remove);
                            
                        end
                        
                        %conversion, in real unit
                        VReal=ADC2Real(BinData,Info.Range(numchan(1,1))/2,Info.Gain(numchan(1,1)),Info.Offset(numchan(1,1)));
                        
                    elseif numchan(1,1)~=numchan(1,2) %if substraction 
                        
                        if ProcessInfo{nprocess,8}==1 %if ICA correction
                            AllBinData=ExtractContinuousData(Info.FilesDir,Info,[],AbsStart,AbsEnd,PerSat2Remove);
                            AllBinDataICA=AllBinData;
                            
                            Mat_ICA_corr=Info.Artef.ICAMat(:,:,2);

                            % %get the chan of interest and chan 2 correct
                            ChanOfInterest=eval(Info.Artef.ChanOfInterest); 
                            ChanToCorrect=eval(Info.Artef.ChanToCorrect);

                            AllBinDataICA(ChanOfInterest,:) = (pinv(Mat_ICA_corr)*(Mat_ICA_corr* AllBinData(ChanOfInterest,:)));

                            %get the value not corrected
                            AllBinData(ChanToCorrect,:)=AllBinDataICA(ChanToCorrect,:);

                            BinData1=AllBinData(numchan(1,1),:);
                            BinData2=AllBinData(numchan(1,2),:);
                            
                            clear AllBinDataICA AllBinData;
                        else
                            BinData1=ExtractContinuousData(Info.FilesDir,Info,numchan(1,1),AbsStart,AbsEnd,PerSat2Remove);
                            BinData2=ExtractContinuousData(Info.FilesDir,Info,numchan(1,2),AbsStart,AbsEnd,PerSat2Remove);
                            
                        end
                        %conversion, in real unit
                        VReal1=ADC2Real(BinData1,Info.Range(numchan(1,1))/2,Info.Gain(numchan(1,1)),Info.Offset(numchan(1,1)));

                        %conversion, in real unit
                        VReal2=ADC2Real(BinData2,Info.Range(numchan(1,2))/2,Info.Gain(numchan(1,2)),Info.Offset(numchan(1,2)));

                        clear BinData1 BinData2;
                        
                        VReal=VReal1-VReal2;

                        clear  VReal1 VReal2;
                    end
                    if size(numchan,1)==2 %if 2 chan for cross corr evaluation 
                         if ProcessInfo{nprocess,8}==1 %if ICA correction
                            AllBinData=ExtractContinuousData(Info.FilesDir,Info,[],AbsStart,AbsEnd,PerSat2Remove);
                            AllBinDataICA=AllBinData;
                            Mat_ICA_corr=Info.Artef.ICAMat(:,:,2);

                            % %get the chan of interest and chan 2 correct
                            ChanOfInterest=eval(Info.Artef.ChanOfInterest); 
                            ChanToCorrect=eval(Info.Artef.ChanToCorrect);

                            AllBinDataICA(ChanOfInterest,:) = (pinv(Mat_ICA_corr)*(Mat_ICA_corr* AllBinData(ChanOfInterest,:)));

                            %get the value not corrected
                            AllBinData(ChanToCorrect,:)=AllBinDataICA(ChanToCorrect,:);
                            %if substraction
                            if numchan(1,1)~=numchan(1,2)
                                BinData1=AllBinData(numchan(1,1),:);
                                BinData1bis=AllBinData(numchan(1,2),:);
                            else
                                BinData1=AllBinData(numchan(1,1),:);
                            end
                            
                            if numchan(2,1)~=numchan(2,2)
                                BinData2=AllBinData(numchan(2,1),:);
                                BinData2bis=AllBinData(numchan(2,2),:);
                            else
                                BinData2=AllBinData(numchan(2,1),:);
                            end
                                
                            
                            
                            clear AllBinDataICA AllBinData;
                         else
                            if numchan(1,1)~=numchan(1,2)%substraction
                                BinData1=ExtractContinuousData(Info.FilesDir,Info,numchan(1,1),AbsStart,AbsEnd,PerSat2Remove);
                                BinData1bis=ExtractContinuousData(Info.FilesDir,Info,numchan(1,2),AbsStart,AbsEnd,PerSat2Remove);
                            else
                                BinData1=ExtractContinuousData(Info.FilesDir,Info,numchan(1,1),AbsStart,AbsEnd,PerSat2Remove);
                            end
                            if numchan(2,1)~=numchan(2,2) %substraction%
                                BinData2=ExtractContinuousData(Info.FilesDir,Info,numchan(2,1),AbsStart,AbsEnd,PerSat2Remove);
                                BinData2bis=ExtractContinuousData(Info.FilesDir,Info,numchan(2,2),AbsStart,AbsEnd,PerSat2Remove);
                            else
                                BinData2=ExtractContinuousData(Info.FilesDir,Info,numchan(2,1),AbsStart,AbsEnd,PerSat2Remove);
                                
                            end
                            
                                
                            
                            
                        end
                        %conversion, in real unit
                        
                        if numchan(1,1)~=numchan(1,2)%substraction
                            VReal1=ADC2Real(BinData1,Info.Range(numchan(1,1))/2,Info.Gain(numchan(1,1)),Info.Offset(numchan(1,1)))- ...
                                ADC2Real(BinData1bis,Info.Range(numchan(1,2))/2,Info.Gain(numchan(1,2)),Info.Offset(numchan(1,2)));
                        else
                            VReal1=ADC2Real(BinData1,Info.Range(numchan(1,1))/2,Info.Gain(numchan(1,1)),Info.Offset(numchan(1,1)));
                            
                        end
                        %conversion, in real unit
                        if numchan(2,1)~=numchan(2,2)%substraction
                            VReal2=ADC2Real(BinData2,Info.Range(numchan(2,1))/2,Info.Gain(numchan(2,1)),Info.Offset(numchan(2,1)))- ...
                                ADC2Real(BinData2bis,Info.Range(numchan(2,2))/2,Info.Gain(numchan(2,2)),Info.Offset(numchan(2,2)));
                        else
                            VReal2=ADC2Real(BinData2,Info.Range(numchan(2,1))/2,Info.Gain(numchan(2,1)),Info.Offset(numchan(2,1)));
                        end

                        %concatenate the 2 chan
                        VReal=[VReal1;VReal2];

                        clear  VReal1 VReal2 BinData BinData2 BinData1bis BinData2bis;
                    end 

                    
                    %replace artefact epoch by NaN
                    


                    %get NaN index and set them to 0
                    idNaN=isnan(VReal);
                    VReal(idNaN)=0; 
                    
                    
                    
                    %check if filtering
                    if isempty(ProcessInfo{nprocess,5})==0 %if filter
                        FilterNum=ProcessInfo{nprocess,5};
                        for nchan=1:size(VReal,1)
                            VReal(nchan,:)=filtfilt(Filt(FilterNum).Hd.sosMatrix,Filt(FilterNum).Hd.ScaleValues,VReal(nchan,:));
                        end
                    end
                    
                    

                    %check if subsampling 
                    if isempty(ProcessInfo{nprocess,4})==0 && ProcessInfo{nprocess,4}>1%if subsampling
                        DecimateFactor=ProcessInfo{nprocess,4};
                        %newFs=Info.Fs/DecimateFactor;
                        for nchan=1:size(VReal,1)
                            VRealTemp(nchan,:) = decimate(VReal(nchan,:),DecimateFactor);
                        end
                        VReal=VRealTemp;
                        clear VRealTemp

                        NewTimeVreal=linspace(AbsStart,AbsEnd,size(VReal,2));
                        newFs=(size(VReal,2)-1)/(AbsEnd-AbsStart);
                        %update idnan 
                        idNaN=idNaN(1:DecimateFactor:end);
                        
                        
                        %if pattern correlation decimate the pattern
                        if ProcessInfo{nprocess,2}==6 %compute the pattern correlation
                   
                                try 
                                    load([Info.FilesDir ProcessInfo{nprocess,3}]);
                                    %convert Pattern to real value
                                    Pattern=ADC2Real(Pattern,Info.Range(numchan(1,1))/2,Info.Gain(numchan(1,1)),Info.Offset(numchan(1,1)));
                                    %use the same decimation as the current chan process    
                                    Pattern= decimate(Pattern,DecimateFactor);

                                catch
                                    Pattern=[];
                                    'error when loading the pattern'
                                end
                        end
                        
                        
                        

                    else
                         
                        NewTimeVreal=linspace(AbsStart,AbsEnd,size(VReal,2));
                        newFs=(size(VReal,2)-1)/(AbsEnd-AbsStart);
                    end

                    %whiteneing the siganl 
                    if length(ProcessInfo(nprocess,:))>=9  && isempty(ProcessInfo{nprocess,9})==0 && ProcessInfo{nprocess,9}==1
                        %add a dc remove (the mean) to avoid problem in
                        %whitening
                        VReal=VReal-mean(VReal(~idNaN));
                        ARModel
                       [VReal,ARModel]=Whitening( VReal,[],1,ARModel,2);
                        
                       
                    end


                    %update NaN into the Vreal variable because decimate and filter couldn't use NaN;
                    VReal(idNaN)=NaN;
%                     
%                     
%                 else %if no need to extract the time lim
%                     if isempty(ProcessInfo{nprocess,4})==0 && ProcessInfo{nprocess,4}>1%if subsampling
%                         DecimateFactor=ProcessInfo{nprocess,4};
%                         NewTimeVreal=linspace(AbsStart,AbsEnd,size(VReal,2));
%                         newFs=size(VReal,2)/(AbsEnd-AbsStart);
%                     else
%                         NewTimeVreal=linspace(AbsStart,AbsEnd,size(VReal,2));
%                         newFs=size(VReal,2)/(AbsEnd-AbsStart);
%                     end
  
                end
                
                
                %generate TimeLim
                



                if ProcessInfo{nprocess,2}==1 %compute TimeSpectrum MMT
 
                    %get nan idx
                    idNaN=isnan(VReal);
                    
                    %DC Remove
                    VReal=VReal-mean(VReal(~idNaN));
                    
                    %replace nan by 0
                    VReal(idNaN)=0;
                    
                    %Process the time frequency analysis
                    MTT(ProcessInfo{nprocess,3}).params.Fs=newFs;
                    [S,t,F]=mtspecgramc(VReal,MTT(ProcessInfo{nprocess,3}).Movingwin,MTT(ProcessInfo{nprocess,3}).params);
                    
                     
                     

                    %remove artefact if necessary
                     if isempty(Artef.ArtefTimeDeb)==0 && isempty(Artef.ArtefTimeEnd)==0 &&  ProcessInfo{nprocess,7}==1
                         for nArtef=1:length(Artef.ArtefTimeDeb)
                             id=find(t+AbsStart>=Artef.ArtefTimeDeb(nArtef) & t+AbsStart<=Artef.ArtefTimeEnd(nArtef));
                             S(id,:)=NaN;
                         end
                     end
                     
                     %replace 0 in TF by Nan
                     if sum(idNaN)>1
                         idTFZero=sum(S,2)==0;
                         S(idTFZero,:)=NaN;
                     end
                         

                    % update the output variable
                     TF.Chan(numchan(1,1),numchan(1,2)).time=t+AbsStart;
                     TF.Chan(numchan(1,1),numchan(1,2)).S=S;
                     TF.Chan(numchan(1,1),numchan(1,2)).F=F;
                     TF.Chan(numchan(1,1),numchan(1,2)).Chnum=ProcessInfo{nprocess,1};
                     TF.Chan(numchan(1,1),numchan(1,2)).params=MTT(ProcessInfo{nprocess,3});
               
                     %plot the Process
                    if FigSave==1
                         h(nprocess+1,nCutPeriod)=subplot('position',[Xo,1-(Yo+(nprocess+1)*(H+Sp)),W,H],'parent',hf(nCutPeriod));
                         
                         if ProcessInfo{nprocess,9}==0 %no whitening
                            pcolor(TF.Chan(numchan(1,1),numchan(1,2)).time./60, TF.Chan(numchan(1,1),numchan(1,2)).F, 10*log10(real(S)'));
                         else % whitening applied
                             pcolor(TF.Chan(numchan(1,1),numchan(1,2)).time./60, TF.Chan(numchan(1,1),numchan(1,2)).F, S');
                         end
                         shading('flat');colormap(jet(2048));
                         ylabel(sprintf('%s', VarName{nprocess}));

                         %set Color setup if necessary
                         if isempty(ProcessInfo{nprocess,6})==0
                            set(h(nprocess+1,nCutPeriod),'Clim',ProcessInfo{nprocess,6});
                         else
                             nprocess
                             VarName{nprocess}
                             get(h(nprocess+1,nCutPeriod),'Clim')
                         end
                    end



                elseif ProcessInfo{nprocess,2}==2 %compute the heart rate
                    %remove artefact if necessary
                    
                     if isempty(Artef.ArtefTimeDeb)==0 && isempty(Artef.ArtefTimeEnd)==0 &&  ProcessInfo{nprocess,7}==1
                         for nArtef=1:length(Artef.ArtefTimeDeb)
                             id=find(NewTimeVreal>=Artef.ArtefTimeDeb(nArtef) & NewTimeVreal<=Artef.ArtefTimeEnd(nArtef));
                             VReal(id)=NaN;
                         end
                     end
                    %Process the heart rate extraction in bpm

                     %find peaks with a min pike high at 1 and maximm bpm is 60 bpm = 60/(300/FS)
                     %warning findpeaks is the matlab function youneed to remove chronux from
                     %the set path
                     [pks,loc] =findpeaks(VReal,'MINPEAKDISTANCE',round(ProcessInfo{nprocess,3}(1)*newFs),'MINPEAKHEIGHT',ProcessInfo{nprocess,3}(2));

                     HB.Chan(numchan(1,1),numchan(1,2)).Peak=pks;
                     HB.Chan(numchan(1,1),numchan(1,2)).timePeak=NewTimeVreal(loc);
                     
                     %compute the bpm
                     HB.Chan(numchan(1,1),numchan(1,2)).BPM=1./diff(NewTimeVreal(loc))*60;
                     
                     %get the time of the instantaneous HB, TBPM is between
                     %each peak detected
                     BPMidx=round((loc(1:end-1)+loc(2:end))/2);
                     
                     
                     %remove Peak detected between file transition (NaN)
                     idBPM2Remove=find(idNaN(BPMidx)==1);
                    
                     BPMidx(idBPM2Remove)=[];
                     HB.Chan(numchan(1,1),numchan(1,2)).BPM(idBPM2Remove)=[];
                     
                     
                     
                     
                     HB.Chan(numchan(1,1),numchan(1,2)).tBPM=NewTimeVreal(BPMidx);
                     
                     
                     %interpolate to have continuous data
                     HBinterp = interp1(HB.Chan(numchan(1,1),numchan(1,2)).tBPM,HB.Chan(numchan(1,1),numchan(1,2)).BPM,NewTimeVreal,'pchip','extrap');
                     
                     
                     %add a low pass at 0.01 hz order 3
                     hp=fdesign.lowpass('N,F3dB',2,0.05,newFs);
                     Hd = design(hp,'butter');
                     
                     
                     %smooth
                     HB.Chan(numchan(1,1),numchan(1,2)).HBsmooth=fastsmooth(HBinterp,newFs*10,3);
                     
                     
                     %HB.Chan(numchan(1,1),numchan(1,2)).HBsmooth = filtfilt(Hd.sosMatrix,Hd.ScaleValues,HBinterp);
                     HB.Chan(numchan(1,1),numchan(1,2)).time=NewTimeVreal;

                     %remove artefact if necessary
                    if isempty(Artef.ArtefTimeDeb)==0 && isempty(Artef.ArtefTimeDeb)==0  &&  ProcessInfo{nprocess,7}==1
                       for nArtef=1:length(Artef.ArtefTimeDeb)
                           id=find(NewTimeVreal>=Artef.ArtefTimeDeb(nArtef) & NewTimeVreal<=Artef.ArtefTimeEnd(nArtef));
                           HB.Chan(numchan(1,1),numchan(1,2)).HBsmooth(id)=NaN;
                       end
                    end

                    %plot the data
                    if FigSave==1
                        h(nprocess+1,nCutPeriod)=subplot('position',[Xo,1-(Yo+(nprocess+1)*(H+Sp)),W,H],'parent',hf(nCutPeriod));
                        plot( HB.Chan(numchan(1,1),numchan(1,2)).tBPM./60, HB.Chan(numchan(1,1),numchan(1,2)).BPM,'b.');hold on;
                        plot( HB.Chan(numchan(1,1),numchan(1,2)).time./60, HB.Chan(numchan(1,1),numchan(1,2)).HBsmooth,'g','LineWidth',2);
                        ylabel(sprintf('%s',VarName{nprocess}));
                        grid on;

                        %set ylim if specify
                        if isempty(ProcessInfo{nprocess,6})==0
                            set(h(nprocess+1,nCutPeriod),'ylim',ProcessInfo{nprocess,6});
                         else
                             nprocess
                             VarName{nprocess}
                             get(h(nprocess+1,nCutPeriod),'ylim')
                         end
                    end



                elseif ProcessInfo{nprocess,2}==3 %raw signal
                    

                    %remove artefact if specify
                    if isempty(Artef.ArtefTimeDeb)==0 && isempty(Artef.ArtefTimeDeb)==0  &&  ProcessInfo{nprocess,7}==1
                         for nArtef=1:length(Artef.ArtefTimeDeb)
                             id=find(NewTimeVreal>=Artef.ArtefTimeDeb(nArtef) & NewTimeVreal<=Artef.ArtefTimeEnd(nArtef));
                             VReal(id)=NaN;
                         end
                    end

                    Raw.Chan(numchan(1,1),numchan(1,2)).Data = VReal;
                    Raw.Chan(numchan(1,1),numchan(1,2)).time =NewTimeVreal;
                    
                    %plot the raw channel  
                    if FigSave==1
                        h(nprocess+1,nCutPeriod)=subplot('position',[Xo,1-(Yo+(nprocess+1)*(H+Sp)),W,H],'parent',hf(nCutPeriod));
                        plot(NewTimeVreal./60,VReal);ylabel(sprintf('%s',VarName{nprocess}));
                        grid on;
                        if isempty(ProcessInfo{nprocess,6})==0
                            set(h(nprocess+1,nCutPeriod),'ylim',ProcessInfo{nprocess,6});
                         else
                             nprocess
                             VarName{nprocess}
                             get(h(nprocess+1,nCutPeriod),'ylim')'
                        end
                    end

                 elseif ProcessInfo{nprocess,2}==4 %peak frequency
                     
                     bandLim=ProcessInfo{nprocess,3};
                     if isfield(TF.Chan(numchan(1,1),numchan(1,2)),'S') && isempty(TF.Chan(numchan(1,1),numchan(1,2)).S)==0
                         t=TF.Chan(numchan(1,1),numchan(1,2)).time;
                         if isempty(Artef.ArtefTimeDeb)==0 && isempty(Artef.ArtefTimeEnd)==0 &&  ProcessInfo{nprocess,7}==1
                             for nArtef=1:length(Artef.ArtefTimeDeb)
                                 id=find(t>=Artef.ArtefTimeDeb(nArtef) & t<=Artef.ArtefTimeEnd(nArtef));
                                 TF.Chan(numchan(1,1),numchan(1,2)).S(id,:)=NaN;
                             end
                         end
                         
                         
                         idband=TF.Chan(numchan(1,1),numchan(1,2)).F>=bandLim(1) & TF.Chan(numchan(1,1),numchan(1,2)).F<bandLim(2);
                         [Smax,iMax]=max(TF.Chan(numchan(1,1),numchan(1,2)).S(:,idband),[],2);
                         
                         %not logical indexes
                         idxband=find(idband==1); 
                         
                         for nbinTF=1:length(TF.Chan(numchan(1,1),numchan(1,2)).time)
                            PeakFqcy(numchan(1,1),numchan(1,2)).Fmax(nbinTF)=TF.Chan(numchan(1,1),numchan(1,2)).F(idxband(iMax(nbinTF)));
                         end
                          PeakFqcy(numchan(1,1),numchan(1,2)).BandLim=bandLim;
                          PeakFqcy(numchan(1,1),numchan(1,2)).time=TF.Chan(numchan(1,1),numchan(1,2)).time;
                          
                           %plot the band power evolution
                        if FigSave==1
                            h(nprocess+1,nCutPeriod)=subplot('position',[Xo,1-(Yo+(nprocess+1)*(H+Sp)),W,H],'parent',hf(nCutPeriod));
                            plot(PeakFqcy(numchan(1,1),numchan(1,2)).time./60,PeakFqcy(numchan(1,1),numchan(1,2)).Fmax);ylabel(sprintf('%s',VarName{nprocess}));hold on;
                            
                           
%                         %add a low pass at 0.01 hz order 3
%                         currFs=1/diff(PeakFqcy(numchan(1,1),numchan(1,2)).time(1:2));
%                          hp=fdesign.lowpass('N,F3dB',2,0.05,currFs);
%                          Hd = design(hp,'butter');
%                          xPeak=filtfilt(Hd.sosMatrix,Hd.ScaleValues,PeakFqcy(numchan(1,1),numchan(1,2)).Fmax);
%                          xPeak2 = vanherk(PeakFqcy(numchan(1,1),numchan(1,2)).Fmax,round(currFs*10),'max','same');
%                         
%                         
%                         %xPeak = hilbert(PeakFqcy(numchan(1,1),numchan(1,2)).Fmax,length(PeakFqcy(numchan(1,1),numchan(1,2)).Fmax));
%                         
%                         
%                         plot(PeakFqcy(numchan(1,1),numchan(1,2)).time./60,PeakFqcy(numchan(1,1),numchan(1,2)).Fmax);ylabel(sprintf('%s',VarName{nprocess}));hold on;
%                         plot(PeakFqcy(numchan(1,1),numchan(1,2)).time./60,xPeak,'g');
%                         plot(PeakFqcy(numchan(1,1),numchan(1,2)).time(1:length(xPeak2))./60,xPeak2,'r');
%                         
                        
                            grid on;
                            if isempty(ProcessInfo{nprocess,6})==0
                                set(h(nprocess+1,nCutPeriod),'ylim',ProcessInfo{nprocess,6});
                             else
                                 nprocess
                                 VarName{nprocess}                                 
                                 set(h(nprocess+1,nCutPeriod),'ylim',bandLim);
                                 get(h(nprocess+1,nCutPeriod),'ylim')'
                            end
                        end
                          
                          
                         
                         
                     end

                     %check if TF exist for this channel
 

                 elseif ProcessInfo{nprocess,2}==5 %band Mean in a band range
                     
                     bandLim=ProcessInfo{nprocess,3};
                     if isfield(TF.Chan(numchan(1,1),numchan(1,2)),'S') && isempty(TF.Chan(numchan(1,1),numchan(1,2)).S)==0
                          t=TF.Chan(numchan(1,1),numchan(1,2)).time;
                         if isempty(Artef.ArtefTimeDeb)==0 && isempty(Artef.ArtefTimeEnd)==0 &&  ProcessInfo{nprocess,7}==1
                             for nArtef=1:length(Artef.ArtefTimeDeb)
                                 id=find(t>=Artef.ArtefTimeDeb(nArtef) & t<=Artef.ArtefTimeEnd(nArtef));
                                 TF.Chan(numchan(1,1),numchan(1,2)).S(id,:)=NaN;
                             end
                         end
                         
                         
                         idband=TF.Chan(numchan(1,1),numchan(1,2)).F>=bandLim(1) & TF.Chan(numchan(1,1),numchan(1,2)).F<bandLim(2);
                         
                         if (size(BandPW,1)>=numchan(1,1) && size(BandPW,2)>=numchan(1,2)) && isempty(BandPW(numchan(1,1),numchan(1,2)))==0
                             nBand=length(BandPW(numchan(1,1),numchan(1,2)).numBand)+1;
                         else
                             nBand=1;
                         end

                         BandPW(numchan(1,1),numchan(1,2)).numBand(nBand).Data=nanmean(TF.Chan(numchan(1,1),numchan(1,2)).S(:,idband),2)';
                         BandPW(numchan(1,1),numchan(1,2)).numBand(nBand).BandLim=bandLim;
                         BandPW(numchan(1,1),numchan(1,2)).numBand(nBand).time=TF.Chan(numchan(1,1),numchan(1,2)).time;
                         
                         
                       %plot the band power evolution
                        if FigSave==1
                            h(nprocess+1,nCutPeriod)=subplot('position',[Xo,1-(Yo+(nprocess+1)*(H+Sp)),W,H],'parent',hf(nCutPeriod));
                            plot(BandPW(numchan(1,1),numchan(1,2)).numBand(nBand).time./60,BandPW(numchan(1,1),numchan(1,2)).numBand(nBand).Data);
                            hold on;  


                            currFs=1/diff(BandPW(numchan(1,1),numchan(1,2)).numBand(nBand).time(1:2));
                             hp=fdesign.lowpass('N,F3dB',2,0.05,currFs);
                             Hd = design(hp,'butter');
                             xPBand=filtfilt(Hd.sosMatrix,Hd.ScaleValues,BandPW(numchan(1,1),numchan(1,2)).numBand(nBand).Data);
                             %plot(BandPW(numchan(1,1),numchan(1,2)).numBand(nBand).time./60,xPBand,'r');
                             ylabel(sprintf('%s',VarName{nprocess}));

                            grid on;
                            if isempty(ProcessInfo{nprocess,6})==0
                                set(h(nprocess+1,nCutPeriod),'ylim',ProcessInfo{nprocess,6});
                             else
                                 nprocess
                                 VarName{nprocess}
                                 get(h(nprocess+1,nCutPeriod),'ylim')'
                            end
                        end

                         
                         
                         
                     end
                     
               elseif ProcessInfo{nprocess,2}==6 %compute the pattern correlation
                   
                   try isempty(Pattern)
                       
                       

                       
                    [corr,l]=xcorr(VReal,Pattern,'none'); 
                    corr(1:length(VReal)-1)=[];
                    
%                     corr2=normxcorr2(Pattern,VReal);
%                     corr2(1:floor(length(Pattern)/2)-1)=[];
%                     corr2(end-floor(length(Pattern)/2)+1:end)=[];
                       
                    %remove artefact if specify
                    if isempty(Artef.ArtefTimeDeb)==0 && isempty(Artef.ArtefTimeDeb)==0  &&  ProcessInfo{nprocess,7}==1
                         for nArtef=1:length(Artef.ArtefTimeDeb)
                             id=find(NewTimeVreal>=Artef.ArtefTimeDeb(nArtef) & NewTimeVreal<=Artef.ArtefTimeEnd(nArtef));
                             corr(id)=NaN;
                         end
                    end 

                    %plot the raw channel 
                    if FigSave==1
                        h(nprocess+1,nCutPeriod)=subplot('position',[Xo,1-(Yo+(nprocess+1)*(H+Sp)),W,H],'parent',hf(nCutPeriod));
                        plot(NewTimeVreal./60,corr);ylabel(sprintf('%s',VarName{nprocess}));
                        grid on;
                        if isempty(ProcessInfo{nprocess,6})==0
                            set(h(nprocess+1,nCutPeriod),'ylim',ProcessInfo{nprocess,6});
                         else
                             nprocess
                             VarName{nprocess}
                             get(h(nprocess+1,nCutPeriod),'ylim')'
                        end
                    end
       
                       
                   catch
                       'Pattern not load'
                   end
                   
                 elseif ProcessInfo{nprocess,2}==7 %ratio of band sum 
                    
                     
                     bandLim=ProcessInfo{nprocess,3};
                     if isfield(TF.Chan(numchan(1,1),numchan(1,2)),'S') && isempty(TF.Chan(numchan(1,1),numchan(1,2)).S)==0
                          t=TF.Chan(numchan(1,1),numchan(1,2)).time;
                         if isempty(Artef.ArtefTimeDeb)==0 && isempty(Artef.ArtefTimeEnd)==0 &&  ProcessInfo{nprocess,7}==1
                             for nArtef=1:length(Artef.ArtefTimeDeb)
                                 id=find(t>=Artef.ArtefTimeDeb(nArtef) & t<=Artef.ArtefTimeEnd(nArtef));
                                 TF.Chan(numchan(1,1),numchan(1,2)).S(id,:)=NaN;
                             end
                         end
                         
                         idband1=TF.Chan(numchan(1,1),numchan(1,2)).F>=bandLim(1) & TF.Chan(numchan(1,1),numchan(1,2)).F<bandLim(2);
                         idband2=TF.Chan(numchan(1,1),numchan(1,2)).F>=bandLim(3) & TF.Chan(numchan(1,1),numchan(1,2)).F<bandLim(4);
                         
                         if (size(RatioBandPW,1)>=numchan(1,1) && size(RatioBandPW,2)>=numchan(1,2)) && isempty(RatioBandPW(numchan(1,1),numchan(1,2)))==0
                             nBand=length(RatioBandPW(numchan(1,1),numchan(1,2)).numBand)+1;
                         else
                             nBand=1;
                         end
% 
%                          RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).Data=nansum(TF.Chan(numchan(1,1),numchan(1,2)).S(:,idband1),2)'./nansum(TF.Chan(numchan(1,1),numchan(1,2)).S(:,idband2),2)';
%                          RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).BandLim=bandLim;
%                          RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).time=TF.Chan(numchan(1,1),numchan(1,2)).time;
%                          
                         
                         SB1=nanmean(TF.Chan(numchan(1,1),numchan(1,2)).S(:,idband1),2)';%.*mean(TF.Chan(numchan(1,1),numchan(1,2)).F(idband1));
                         SB2=nanmean(TF.Chan(numchan(1,1),numchan(1,2)).S(:,idband2),2)';%.*mean(TF.Chan(numchan(1,1),numchan(1,2)).F(idband2));
%                          
%                          SB1=nanmax(TF.Chan(numchan(1,1),numchan(1,2)).S(:,idband1),[],2)';%.*mean(TF.Chan(numchan(1,1),numchan(1,2)).F(idband1));
%                          SB2=nanmax(TF.Chan(numchan(1,1),numchan(1,2)).S(:,idband2),[],2)';%.*mean(TF.Chan(numchan(1,1),numchan(1,2)).F(idband2));
% %                          

                         %SB1=nanmean(TF.Chan(numchan(1,1),numchan(1,2)).S(:,idband1).*repmat(TF.Chan(numchan(1,1),numchan(1,2)).F(idband1),size(TF.Chan(numchan(1,1),numchan(1,2)).S,1),1),2);
                         %SB2=nanmean(TF.Chan(numchan(1,1),numchan(1,2)).S(:,idband2).*repmat(TF.Chan(numchan(1,1),numchan(1,2)).F(idband2),size(TF.Chan(numchan(1,1),numchan(1,2)).S,1),1),2);

                         RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).Data=SB1./SB2;
                         RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).BandLim=bandLim;
                         RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).time=TF.Chan(numchan(1,1),numchan(1,2)).time;
                         

                         
                         %replace NaN per Zeros
                         idNan=isnan(RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).Data);
                         RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).Data(idNan)=0;
                         
                       
                        
                        
                        currFs=1/diff(RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).time(1:2));
                         hp=fdesign.lowpass('N,F3dB',2,0.05,currFs);
                         Hd = design(hp,'butter');
                         xPBand=filtfilt(Hd.sosMatrix,Hd.ScaleValues,RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).Data);

                        %plot the band power evolution
                        if FigSave==1
                            h(nprocess+1,nCutPeriod)=subplot('position',[Xo,1-(Yo+(nprocess+1)*(H+Sp)),W,H],'parent',hf(nCutPeriod));
                            plot(RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).time./60,RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).Data);
                            hold on;
                            % plot(RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).time./60,xPBand,'r');
                             ylabel(sprintf('%s',VarName{nprocess}));








                            grid on;
                            if isempty(ProcessInfo{nprocess,6})==0
                                set(h(nprocess+1,nCutPeriod),'ylim',ProcessInfo{nprocess,6});
                             else
                                 nprocess
                                 VarName{nprocess}
                                 get(h(nprocess+1,nCutPeriod),'ylim')'
                            end
                    end

                         
                         
                         
                     end
                     
                 elseif ProcessInfo{nprocess,2}==8 %peak detection
                     
                     %remove artefact if necessary

                     if isempty(Artef.ArtefTimeDeb)==0 && isempty(Artef.ArtefTimeEnd)==0 &&  ProcessInfo{nprocess,7}==1
                         for nArtef=1:length(Artef.ArtefTimeDeb)
                             id=find(NewTimeVreal>=Artef.ArtefTimeDeb(nArtef) & NewTimeVreal<=Artef.ArtefTimeEnd(nArtef));
                             VReal(id)=NaN;
                         end
                     end

                     %find peaks with a min pike high
                     %warning findpeaks is the matlab function youneed to remove chronux from
                     
                     MeanVreal=nanmean(VReal);
                     StdVreal=nanstd(VReal);
                     NbStd=ProcessInfo{nprocess,3}(2);

                    
                    [pks,loc] =findpeaks(sign(ProcessInfo{nprocess,3}(2))*VReal,'MINPEAKDISTANCE',round(ProcessInfo{nprocess,3}(1)*newFs),'MINPEAKHEIGHT',sign(ProcessInfo{nprocess,3}(2))*MeanVreal+abs(NbStd)*StdVreal);

                    %plot(NewTimeVreal./60,VReal);
                    if FigSave==1
                        h(nprocess+1,nCutPeriod)=subplot('position',[Xo,1-(Yo+(nprocess+1)*(H+Sp)),W,H],'parent',hf(nCutPeriod));
                        ylabel(sprintf('%s',VarName{nprocess}));
                        grid on;
                        hold on;

                        plot([NewTimeVreal(1)./60 NewTimeVreal(end)./60],[MeanVreal+NbStd*StdVreal MeanVreal+NbStd*StdVreal],'g')
                        hold on;
                        plot(NewTimeVreal(loc)./60,VReal(loc),'k*');

                        Peak(numchan(1,1),numchan(1,2)).PeakVal=sign(ProcessInfo{nprocess,3}(2))*pks;
                        Peak(numchan(1,1),numchan(1,2)).PeakTime=NewTimeVreal(loc);
                        Peak(numchan(1,1),numchan(1,2)).MinPeakDist=round(ProcessInfo{nprocess,3}(1)*newFs);
                        Peak(numchan(1,1),numchan(1,2)).MinPeakThresh=MeanVreal+NbStd*StdVreal;
                        Peak(numchan(1,1),numchan(1,2)).MeanVal=MeanVreal;
                        Peak(numchan(1,1),numchan(1,2)).StdVal=StdVreal;
                        Peak(numchan(1,1),numchan(1,2)).NbStd=NbStd;

                        if isempty(ProcessInfo{nprocess,6})==0
                            set(h(nprocess+1,nCutPeriod),'ylim',ProcessInfo{nprocess,6});
                        else
                             nprocess
                             VarName{nprocess}
                             get(h(nprocess+1,nCutPeriod),'ylim')'
                        end
                    end
                    
                    
                    
                 elseif ProcessInfo{nprocess,2}==9 %XCorr and peak detection from covariance
                     
                     %remove artefact if necessary

                    
                    if Pass==1
                        
                         if isempty(Artef.ArtefTimeDeb)==0 && isempty(Artef.ArtefTimeEnd)==0 &&  ProcessInfo{nprocess,7}==1
                             for nArtef=1:length(Artef.ArtefTimeDeb)
                                 id=find(NewTimeVreal>=Artef.ArtefTimeDeb(nArtef) & NewTimeVreal<=Artef.ArtefTimeEnd(nArtef));
                                 VReal(:,id)=NaN;
                             end
                         end
                        
                        
                        
                        Win=ProcessInfo{nprocess,3}(1);
                        WinPt=round(Win*newFs);
                        Step=ProcessInfo{nprocess,3}(2);
                        StepPt=newFs*Step; %step in sample
                        XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Win=Win;
                        XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Step=Step;
                        
                        MaxVal=size(VReal,2)-mod(size(VReal,2),WinPt)-WinPt;
                        Currstep=0;
                        barreXCor=waitbar(0,'Cross corr evaluation in progress');
                        %generate the step list
                        
                        
                                       
                        [lag_time,twin,xcl]=timewindow_xcorr(VReal(1,:),VReal(2,:),newFs, Win,Step,0,1);
                        
                        XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time=linspace(AbsStart+Win/2,AbsEnd-Win/2,length(xcl));
                        XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Data=xcl';
                        XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).NewFs=newFs;

                        
                        
%                         Steplist=round([1:StepPt:MaxVal]);
%                         for nn=1:length(Steplist)
%                            n=Steplist(nn);
%                            Xpart=VReal(1,n:n+WinPt-1);
%                            Ypart=VReal(2,n:n+WinPt-1);
%                            pXCor=cov(Xpart',Ypart'); 
%                            Currstep=Currstep+1;
%                            XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time(Currstep)=NewTimeVreal(round(n+(Step/2)));
%                            XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Data(Currstep)=pXCor(2,1);
%                            waitbar(n/MaxVal,barreXCor);
%                         end
                        close(barreXCor)
                        
                        %plot data
                         if FigSave==1
                            h(nprocess+1,nCutPeriod)=subplot('position',[Xo,1-(Yo+(nprocess+1)*(H+Sp)),W,H],'parent',hf(nCutPeriod));

                            ylabel(sprintf('%s',VarName{nprocess}));


                            plot(XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time/60,XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Data,'r');
                            grid on;
                            if isempty(ProcessInfo{nprocess,6})==0
                                set(h(nprocess+1,nCutPeriod),'ylim',ProcessInfo{nprocess,6});
                            else
                                 nprocess
                                 VarName{nprocess}
                                 get(h(nprocess+1,nCutPeriod),'ylim')'
                            end 
                        end
                     
                    elseif Pass==2
                        NewTimeVreal=XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time;
                        newFs=(length(NewTimeVreal)-1)/(NewTimeVreal(end)-NewTimeVreal(1));
                        if isempty(Artef.ArtefTimeDeb)==0 && isempty(Artef.ArtefTimeEnd)==0 &&  ProcessInfo{nprocess,7}==1
                             for nArtef=1:length(Artef.ArtefTimeDeb)
                                 id=find(NewTimeVreal>=Artef.ArtefTimeDeb(nArtef) & NewTimeVreal<=Artef.ArtefTimeEnd(nArtef));
                                 XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Data(id)=NaN;
                             end
                         end
                   
                         %find peaks with a min pike high
                         %warning findpeaks is the matlab function youneed to remove chronux from
                         MinPeakDist=ProcessInfo{nprocess,3}(3);
                         MinPeakThreshold=ProcessInfo{nprocess,3}(4);
                         MeanXCor=nanmean(XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Data);
                         StdXCor=nanstd(XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Data);
                         NbStd=MinPeakThreshold;

                         
%                          HalfPeakWidth=round(MinPeakDist*newFs/2);
%                          smoothtype=1;
%                          SlopeThreshold=0.7*(HalfPeakWidth*2)^-2;
%                          AmpThreshold=abs(NbStd)*StdXCor;
%                          smoothwidth=round(HalfPeakWidth);
%                          FitWidth=smoothwidth;

                         HalfPeakWidth=round(MinPeakDist*newFs/2);
%                          smoothtype=1;
%                          SlopeThreshold=0.7*(HalfPeakWidth*2)^-2;
                          AmpThreshold=abs(NbStd)*StdXCor+MeanXCor;
%                          smoothwidth=round(HalfPeakWidth);
%                          FitWidth=smoothwidth*2;
                         
                         %inverse the XCorr if look for negative correlation
                         YNaN=sign(MinPeakThreshold)*XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Data;
                         X=1:length(YNaN);
                         
                         %RemoveNaN
                         IdNotNaN=isnan(YNaN)==0;
                         if sum(IdNotNaN)>=1
                         %interp to fill the Blank with the nearestneighbor
                            Y = interp1(X(IdNotNaN),YNaN(IdNotNaN),X,'nearest','extrap');
                         else
                             Y=YNaN;
                         end
                         
                         %remove the negative correlation
                        %Y(Y<=0)=0;
                        %SlopeThreshold,AmpThreshold,smoothwidth,FitWidth,smoothtype
                        %P=findpeaksG(X,Y,SlopeThreshold,AmpThreshold,smoothwidth,FitWidth,smoothtype);
                        [P,loc]=findpeaks(Y,'minpeakheight',AmpThreshold,'minpeakdistance',HalfPeakWidth.*2);
                        
%                         if sum(P(1,:))==0
%                             loc=[];
%                         else
%                             loc=round(P(:,2));
%                         end
                        
                        %%remove the 
                        %inverse the Peak value if looking for negative XCorr
                        pks=sign(MinPeakThreshold)*Y(loc);
                       %[pks,loc] =findpeaks(sign(MinPeakThreshold)*XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Data,'MINPEAKDISTANCE',round(MinPeakDist/XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Step),'MINPEAKHEIGHT',sign(MinPeakThreshold)*MeanXCor+abs(NbStd)*StdXCor);

                        %plot(NewTimeVreal./60,VReal);

                        XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).PeakVal=pks;
                        XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).PeakTime=XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time(loc);
                        XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Id=loc;
                        XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).MinPeakDist=MinPeakDist;
                        XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).MinPeakThresh=MinPeakThreshold;
                        XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).MinPeakDistPt=round(MinPeakDist/XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Step);
                        XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).MinPeakThreshVal=MeanXCor+NbStd*StdXCor;
                        XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).MeanVal=MeanXCor;
                        XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).StdVal=StdXCor;
                        XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).NbStd=NbStd;





                        if FigSave==1
                            h(nprocess+1,nCutPeriod)=subplot('position',[Xo,1-(Yo+(nprocess+1)*(H+Sp)),W,H],'parent',hf(nCutPeriod));
                            ylabel(sprintf('%s',VarName{nprocess}));


                            plot(Xnumchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2).time/60,XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Data,'b');hold on;
                            plot(XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time(X)/60,sign(MinPeakThreshold)*Y,'r');
                            plot([XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time(1)./60 XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time(end)./60],[MeanXCor+NbStd*StdXCor MeanXCor+NbStd*StdXCor],'g')
                            grid on;
                           % plot(XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time(loc)./60,XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Data(loc),'k*');
                            plot(XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time(X(loc))/60,sign(MinPeakThreshold)*Y(loc),'k*');

% 
%                         plot(XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time/60,Y,'r');
%                         hold on;
%                         plot([XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time(1)./60 XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time(end)./60],[MeanXCor+abs(NbStd)*StdXCor MeanXCor+abs(NbStd)*StdXCor],'g')
%                         grid on;
%                         plot(XCor(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time(loc)./60,Y(loc),'k*');


    
                            if isempty(ProcessInfo{nprocess,6})==0
                                set(h(nprocess+1,nCutPeriod),'ylim',ProcessInfo{nprocess,6});
                            else
                                 nprocess
                                 VarName{nprocess}
                                 get(h(nprocess+1,nCutPeriod),'ylim')'
                            end  
                        end
                    end
                        
                        
                elseif ProcessInfo{nprocess,2}==10 %hilbert real part
                   

                    
                    
                    
                    %filtering
                    % design the filter
                    bandLim=ProcessInfo{nprocess,3};
                    HilbFilt.hp  = fdesign.bandpass('N,Fc1,Fc2',10,bandLim(1),bandLim(2),newFs);
                    HilbFilt.Hd = design(HilbFilt.hp,'butter');
                    
                    VrealFilt=filtfilt(HilbFilt.Hd.sosMatrix,HilbFilt.Hd.ScaleValues,VReal);
                    HilbData=real(hilbert(VrealFilt));
                    %remove artefact if specify
                    if isempty(Artef.ArtefTimeDeb)==0 && isempty(Artef.ArtefTimeDeb)==0  &&  ProcessInfo{nprocess,7}==1
                         for nArtef=1:length(Artef.ArtefTimeDeb)
                             id=find(NewTimeVreal>=Artef.ArtefTimeDeb(nArtef) & NewTimeVreal<=Artef.ArtefTimeEnd(nArtef));
                             HilbData(id)=NaN;
                         end
                    end
                    
                    HilbertReal(numchan(1,1),numchan(1,2)).numBand(1).Data=HilbData;
                    HilbertReal(numchan(1,1),numchan(1,2)).numBand(1).BandLim=bandLim;
                    HilbertReal(numchan(1,1),numchan(1,2)).numBand(1).time=NewTimeVreal;
                    
                    

                    %plot the raw channel
                     if FigSave==1
                         h(nprocess+1,nCutPeriod)=subplot('position',[Xo,1-(Yo+(nprocess+1)*(H+Sp)),W,H],'parent',hf(nCutPeriod));
                        plot(NewTimeVreal./60,HilbData);ylabel(sprintf('%s',VarName{nprocess}));
                        grid on;
                        if isempty(ProcessInfo{nprocess,6})==0
                            set(h(nprocess+1,nCutPeriod),'ylim',ProcessInfo{nprocess,6});
                         else
                             nprocess
                             VarName{nprocess}
                             get(h(nprocess+1,nCutPeriod),'ylim')'
                        end
                    end
             elseif ProcessInfo{nprocess,2}==11 %ratio of mean PW b1 f1/mean PW(f2+f3)/2
                     
                     bandLim=ProcessInfo{nprocess,3};
                     if isfield(TF.Chan(numchan(1,1),numchan(1,2)),'S') && isempty(TF.Chan(numchan(1,1),numchan(1,2)).S)==0
                          t=TF.Chan(numchan(1,1),numchan(1,2)).time;
                         if isempty(Artef.ArtefTimeDeb)==0 && isempty(Artef.ArtefTimeEnd)==0 &&  ProcessInfo{nprocess,7}==1
                             for nArtef=1:length(Artef.ArtefTimeDeb)
                                 id=find(t>=Artef.ArtefTimeDeb(nArtef) & t<=Artef.ArtefTimeEnd(nArtef));
                                 TF.Chan(numchan(1,1),numchan(1,2)).S(id,:)=NaN;
                             end
                         end
                         
                         idband1=TF.Chan(numchan(1,1),numchan(1,2)).F>=bandLim(1) & TF.Chan(numchan(1,1),numchan(1,2)).F<bandLim(2);
                         idband2=TF.Chan(numchan(1,1),numchan(1,2)).F>=bandLim(3) & TF.Chan(numchan(1,1),numchan(1,2)).F<bandLim(4);
                         idband3=TF.Chan(numchan(1,1),numchan(1,2)).F>=bandLim(5) & TF.Chan(numchan(1,1),numchan(1,2)).F<bandLim(6);
                         
                         if (size(RatioBandPW,1)>=numchan(1,1) && size(RatioBandPW,2)>=numchan(1,2)) && isempty(RatioBandPW(numchan(1,1),numchan(1,2)))==0
                             nBand=length(RatioBandPW(numchan(1,1),numchan(1,2)).numBand)+1;
                         else
                             nBand=1;
                         end
                         
                         SB1=nanmean(TF.Chan(numchan(1,1),numchan(1,2)).S(:,idband1),2)';%.*mean(TF.Chan(numchan(1,1),numchan(1,2)).F(idband1));
                         SB2=nanmean(TF.Chan(numchan(1,1),numchan(1,2)).S(:,idband2|idband3),2)';%.*mean(TF.Chan(numchan(1,1),numchan(1,2)).F(idband2));
                         %SB3=nanmean(TF.Chan(numchan(1,1),numchan(1,2)).S(:,idband3),2)';%.*mean(TF.Chan(numchan(1,1),numchan(1,2)).F(idband3));

                         RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).Data=SB1./SB2;
                         RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).BandLim=bandLim;
                         RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).time=TF.Chan(numchan(1,1),numchan(1,2)).time;
                         
                         %replace NaN per Zeros
                         idNan=isnan(RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).Data);
                         RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).Data(idNan)=0;
                         
                       %plot the band power evolution
                        if FigSave==1
                            h(nprocess+1,nCutPeriod)=subplot('position',[Xo,1-(Yo+(nprocess+1)*(H+Sp)),W,H],'parent',hf(nCutPeriod));
                            plot(RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).time./60,RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).Data);
                            hold on;


                            currFs=1/diff(RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).time(1:2));
                             hp=fdesign.lowpass('N,F3dB',2,0.05,currFs);
                             Hd = design(hp,'butter');
                             xPBand=filtfilt(Hd.sosMatrix,Hd.ScaleValues,RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).Data);
                             %plot(RatioBandPW(numchan(1,1),numchan(1,2)).numBand(nBand).time./60,xPBand,'r');
                             ylabel(sprintf('%s',VarName{nprocess}));








                            grid on;
                            if isempty(ProcessInfo{nprocess,6})==0
                                set(h(nprocess+1,nCutPeriod),'ylim',ProcessInfo{nprocess,6});
                             else
                                 nprocess
                                 VarName{nprocess}
                                 get(h(nprocess+1,nCutPeriod),'ylim')'
                            end
                        end

                         
                         
                         
                     end
              elseif ProcessInfo{nprocess,2}==12 %Coherence
                     
                     %remove artefact if necessary

                    
                    if Pass==1
                        %DC Remove
                        VReal(1,:)=VReal(1,:)-mean(VReal(1,:));
                        VReal(2,:)=VReal(2,:)-mean(VReal(2,:));
                        MTT(ProcessInfo{nprocess,3}).params.Fs=newFs;
                           
                        [C,phi,S12,S1,S2,t,f]=cohgramc(VReal(1,:)',VReal(2,:)',MTT(ProcessInfo{nprocess,3}).Movingwin,MTT(ProcessInfo{nprocess,3}).params);

                                             % update the output variable
                         Cohe.Chan(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time=t+AbsStart;
                         Cohe.Chan(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).C=C;
                         Cohe.Chan(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).F=f;
                         Cohe.Chan(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).phi=phi;
%                          Cohe.Chan(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Cerr=Cerr;
                         
                         Cohe.Chan(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).Chnum=ProcessInfo{nprocess,1};
                         Cohe.Chan(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).params=MTT(ProcessInfo{nprocess,3});
                        

                       %remove artefact if necessary
                         if isempty(Artef.ArtefTimeDeb)==0 && isempty(Artef.ArtefTimeEnd)==0 &&  ProcessInfo{nprocess,7}==1
                             for nArtef=1:length(Artef.ArtefTimeDeb)
                                 id=find(t+AbsStart>=Artef.ArtefTimeDeb(nArtef) & t+AbsStart<=Artef.ArtefTimeEnd(nArtef));
                                 S(id,:)=NaN;
                             end
                         end

 
               
                     %plot the Process
                        if FigSave==1
                             h(nprocess+1,nCutPeriod)=subplot('position',[Xo,1-(Yo+(nprocess+1)*(H+Sp)),W,H],'parent',hf(nCutPeriod));
                             pcolor(Cohe.Chan(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).time./60, Cohe.Chan(numchan(1,1),numchan(1,2),numchan(2,1),numchan(2,2)).F, real(C)');

                             shading('flat');colormap(jet(2048));
                             ylabel(sprintf('%s', VarName{nprocess}));

                             %set Color setup if necessary
                             if isempty(ProcessInfo{nprocess,6})==0
                                set(h(nprocess+1,nCutPeriod),'Clim',ProcessInfo{nprocess,6});
                             else
                                 nprocess
                                 VarName{nprocess}
                                 get(h(nprocess+1,nCutPeriod),'Clim')
                             end
                        end

               
  
                    end

                end
                
                
                
                clear VReal;

                %link all axes and specifiy the x label
                if FigSave==1
                    linkaxes(h(:,nCutPeriod),'x');zoom on;
                    
                    
                     



                    set(h(1,nCutPeriod),'xlim',[AbsStart./60 AbsEnd./60],'FontSize',6);
                    set(h(nprocess+1,nCutPeriod),'xgrid','on','ygrid','on');
                    
                    if nprocess<size(ProcessInfo,1)
                        set(h(nprocess+1,nCutPeriod),'xticklabel',[]);
                     else
                         xtick=get(h(nprocess+1,nCutPeriod),'xtick');
                         
                         if isempty(LightOn)==0
                         
                             T0Date=TimeStartExp+xtick/60/24;

                             %find the first round hour
                             newtick=[];
                             RoundStep=24;

                             %etermine the spce step
                             difftime=0;
                             while difftime==0
                                 EndTime=floor(T0Date(end)*RoundStep)/RoundStep;
                                 DebTime=ceil(T0Date(1)*RoundStep)/RoundStep;
                                 difftime=EndTime-DebTime;
                                 RoundStep=RoundStep/2;
                             end

                             Steptick=1/24;
                             newtick=DebTime:Steptick:EndTime;
                             while length(newtick)>8          
                                newtick=DebTime:Steptick:EndTime; 
                                Steptick=Steptick*2;
                             end
                             newtickstr=datestr(newtick,'dd/mm HH:MM:SS');

                             %define 10 tick
                             xticklabel=datestr(T0Date,'hh:mm:ss');
                             for nax=1:length(h)
                                 set(h(nax,nCutPeriod),'xtick',(newtick-TimeStartExp)*60*24,'xticklabel',newtickstr)
                             end
                         end
                         
                     end
                end



%                 if nprocess~=size(ProcessInfo,1) %add the xtick and xlabel
%                     set(h(nprocess+1,nCutPeriod),'Xticklabel',[]);
%                  else
%                      xlabel('Abs Time');
%                       Xtick=get(h(nprocess+1,nCutPeriod),'xtick');
%                       XTicknew=Xtick/24/60+Info.BinFiles(1).TStart;
%                       formatOut = 'HH:MM:SS.FFF';                     
%                       XTickLabelnew={datestr(XTicknew,formatOut)};
%                       set(h(nprocess+1,nCutPeriod),'xticklabel',XTickLabelnew);
% 
%                 end

                
            catch ME1
                 ME1.identifier
                ME1.message
                for i=1:length(ME1.stack)
                    ME1.stack(i)
                end
            end


        end
        if exist('MainOutDir')==0 || (exist('MainOutDir')~=0 && isempty(MainOutDir)==1)
            CurrDir=Info.FilesDir;
        else
            CurrDir=MainOutDir;
        end
        %save figure
        if FigSave==1
            
            if isempty(LightOn)==0 %if light cycle should be display
                
                AxesLight(nCutPeriod)=subplot('position',[Xo,1-Yo,W,Yo-2*Sp]);
                set(AxesLight(nCutPeriod),'xgrid','off','ygrid','off','xtick',[],'ytick',[]);
                DispDurationMin=AbsEnd./60-AbsStart./60;
                 mintime=floor(T0Date(1)*24)/24; %arrondit inférieur en heure
                 maxtime= ceil(T0Date(1)*24+DispDurationMin/60)/24;%arrondit supérieur en heure
                 VecNum=mintime:1/24:maxtime;
                 a=datevec(VecNum);
                 id=a(:,4)==LightOn(1)| a(:,4)==LightOn(2);
                 Light=a(:,4)>=LightOn(1) & a(:,4)<=LightOn(2)-1;
                 
                 [Xpatch,Ypatch]=GenarateXYPatch(1,[0 length(VecNum)]);
                 LightColor=[0 0 0;1 1 0];

                 patch('XData',Xpatch*60+(mintime-TimeStartExp)*24*60,'YData',Ypatch,'CDataMapping','direct','parent',AxesLight(nCutPeriod),'edgecolor','none','FaceColor','flat','FaceVertexCData',LightColor(Light+1,:));
    
                linkaxes([h(:,nCutPeriod);AxesLight(nCutPeriod)],'x');
                
                
            end
            
            
            
            
            if Pass==1 || isempty(Pass)==1
                TifFileName=sprintf('%s\\Computed Data\\%s %s process %06.1f min-%06.1f min Pass1.tif',CurrDir,OutNamePref,Info.ExpFileName(1:end-4),AbsStart./60,AbsEnd./60);
                FigFileName=sprintf('%s\\Computed Data\\%s %s process %06.1f min-%06.1f minv Pass1.fig',CurrDir,OutNamePref,Info.ExpFileName(1:end-4),AbsStart./60,AbsEnd./60);
            elseif Pass==2
                TifFileName=sprintf('%s\\Computed Data\\%s %s process %06.1f min-%06.1f min Pass2.tif',CurrDir,OutNamePref,Info.ExpFileName(1:end-4),AbsStart./60,AbsEnd./60);
                FigFileName=sprintf('%s\\Computed Data\\%s %s process %06.1f min-%06.1f min Pass2.fig',CurrDir,OutNamePref,Info.ExpFileName(1:end-4),AbsStart./60,AbsEnd./60);
            end
            
            hgsave(hf(nCutPeriod), FigFileName, '-v7.3')
    
            set(hf(nCutPeriod),'PaperPositionMode','Auto');
            print(hf(nCutPeriod),TifFileName,'-dtiff','-r100');
            %saveas(hf(nbin,nCutPeriod),FigFileName);
           close(hf(nCutPeriod));
            
        end


        Header.VarName=VarName;
  


        if MatSave==1
            if Pass==1 || isempty(Pass)==1
                MatFileName=sprintf('%s\\Computed Data\\%s %s process %06.1f min-%06.1f min Pass1.mat',CurrDir,OutNamePref,Info.ExpFileName(1:end-4),T0Cut(nCutPeriod),T0Cut(nCutPeriod)+TotDurationCut(nCutPeriod));
                save(MatFileName,'Header','HB','Raw','TF','Artef','Filt','BandPW','PeakFqcy','RatioBandPW','Peak','XCor','HilbertReal','Cohe','ARModel','-v7.3');
            elseif Pass==2
                MatFileName=sprintf('%s\\Computed Data\\%s %s process %06.1f min-%06.1f min Pass2.mat',CurrDir,OutNamePref,Info.ExpFileName(1:end-4),T0Cut(nCutPeriod),T0Cut(nCutPeriod)+TotDurationCut(nCutPeriod));
                save(MatFileName,'Header','HB','Raw','Artef','Filt','BandPW','PeakFqcy','RatioBandPW','Peak','XCor','HilbertReal','Cohe','ARModel','-v7.3');
            end
            
        end
        
        waitbar(nCutPeriod/length(T0Cut),barre,sprintf('Keep Cool, %3.1f%% processed',nCutPeriod/length(T0Cut)*100));
        
        CurrDir=Info.FilesDir;  
    end
          close(barre);
end




% function SaveData2Mat(Xout,Label,allFileWithoutExt,State2Process,ProcessInfo)
% 
%     %concatenatelabel and Xout
% 
%     Data2Export=num2cell(Xout);
%     Data2Export=cat(1,Label,Data2Export);
% 
%     %save data
%     xlsname=sprintf('%s_Compute.xls',allFileWithoutExt);
%     xlswrite(xlsname,Data2Export);
% 
% 
%     matfile=sprintf('%s_Compute.mat',allFileWithoutExt);
%     save(matfile,'Data2Export','ProcessInfo');
%    
%     warndlg('data saved in a mat file');
% 
% end
% 
% function SaveData2XLS(Xout,Label,allFileWithoutExt,State2Process)
% 
%     %concatenatelabel and Xout
% 
%     Data2Export=num2cell(Xout);
%     Data2Export=cat(1,Label,Data2Export);
% 
%     %save data
%     xlsname=sprintf('%s_Compute.xls',allFileWithoutExt);
%     xlswrite(xlsname,Data2Export);
% 
%     %for each state
% 
%     %Msave in a specific sheet
% 
% 
%     for nstate=1:length(State2Process)
% 
%         XcuurState=[];
%         DataCurrState2Export=[];
%         idxstate=Xout(:,1)==State2Process(nstate);
%         XcurrState=Xout(idxstate,:);
%         DataCurrState2Export=num2cell(XcurrState);
%         DataCurrState2Export=cat(1,Label,DataCurrState2Export);
%         xlswrite(xlsname,DataCurrState2Export,sprintf('State %d',State2Process(nstate)));
%     end
% 
%    
%     warndlg('data saved in an xls file');
% 
% end

% %% load exp file
% function Info=loadEXP
% 
%     try
%         load TempDir
%         [FileName,DefaultDir,FilterIndex]=uigetfile({'*.exp'},'Pick files','multiselect','off',DefaultDir);
%     catch
%         [FileName,DefaultDir,FilterIndex]=uigetfile({'*.exp'},'Pick files','multiselect','off','c:\');
%      end
% 
%      if FilterIndex==1
% 
%         try
% 
%             save TempDir DefaultDir;
%             %save the parameters into Info
% 
%             Info.FilesDir=DefaultDir;
%             Info.ExpFileName=FileName;
%             %create variable with all filename
% 
%             [FileDir,FileWithoutExt,ext] = fileparts(FileName);
%             allFileName = fullfile(DefaultDir,[FileWithoutExt ext]);
%             allFileWithoutExt = fullfile(DefaultDir,FileWithoutExt);
% 
%             %Read all exp file
%             [ s ] = xml2struct(allFileName);
%             Channels=s.Animal.Acquisition.Channels.Channel;
%             ChLabel={};
%             Gain=[];
%             Offset=[];
%             Range=[];
% 
% 
%             for nchan=1:length(Channels)
%                 ChLabel=[ChLabel Channels{nchan}.Name.Text];
%                 Offset=[Offset str2double(Channels{nchan}.Offset.Text)];
%                 Gain=[Gain str2double(Channels{nchan}.Gain.Text)];
%                 Range=[Range str2double(Channels{nchan}.AcquisitionRangeMax.Text).*2];
%             end
% 
% 
%             % Paramètres d'acquisition
%             Info.NbRecChan = str2double(s.Animal.Acquisition.NbChan.Text);
%             Info.ChLabel = ChLabel;
%             Info.Fs = str2double(s.Animal.Acquisition.SamplingRate.Text);
%             Info.Gain = Gain;
%             Info.Range = Range;
%             Info.Offset =Offset; % le point 0 du 16 bit
%             Info.Filename = [FileWithoutExt ext];
% 
%             try %try to importe the states value and color
%                 AnimalMatFileName=[allFileName(1:end-3) 'mat'];
%                 load(AnimalMatFileName);
%                 Info.State=Animal.PStates;
%                  
%                 try %try to load ICA detection 
%                     Info.Artef=[];
%                     Info.Artef=PArtef;
%  
%                 catch
%                     Info.Artef=[];
%                 end
% 
%             catch
%                 Info.State=[];
%                 Info.Artef=[];
%             end
%            
%                 
% 
%             %bin files
%             BinFiles=s.Animal.Acquisition.Files.File;
% 
%             if length(BinFiles)==1 %if only one bin file there is no binfile idx is not a cell
%                 BinFiles={BinFiles}; 
%             end
%             
%             for nbin=1:length(BinFiles)
%                 Info.BinFiles(nbin).FileName=BinFiles{nbin}.FileName.Text;
%                 formatOut = 'yyyy-mm-ddTHH:MM:SS.FFF';
%                 try
%                     Info.BinFiles(nbin).TStart=datenum(BinFiles{nbin}.TStart.Text,formatOut);%Start in relative time
%                 catch
%                     Info.BinFiles(nbin).TStart=datenum([BinFiles{nbin}.TStart.Text '.000'],formatOut);%Start in relative time
%                 end
% 
%                 Info.BinFiles(nbin).Duration=GetBinDuration(fullfile(Info.FilesDir,Info.BinFiles(nbin).FileName),Info.NbRecChan)/Info.Fs;
%             end
% 
%             %Hypno files
%             HypnoFiles=s.Animal.Hypnogram.Files.File;
%             if length(HypnoFiles)==1 %if only one bin file there is no binfile idx is not a cell
%                 HypnoFiles={HypnoFiles};
%             end
% 
% 
%             for nhyp=1:length(HypnoFiles)
%                 Info.HypnoFiles(nhyp).FileName=HypnoFiles{nhyp}.FileName.Text;
%                 formatOut = 'yyyy-mm-ddTHH:MM:SS';
%                 Info.HypnoFiles(nhyp).TStart=datenum(HypnoFiles{nhyp}.TStart.Text,formatOut);%Start in relative time
%                 NumberOfSamples=GetBinDuration( fullfile(Info.FilesDir,Info.HypnoFiles(nhyp).FileName),1);
%                 Info.HypnoFiles(nhyp).Duration=NumberOfSamples;
%             end
% 
% 
%             %reorder the file if not in alphabetical order 
% 
%             [ListName,idx]=sort({Info.BinFiles(:).FileName});
% 
%             BinName={Info.BinFiles(:).FileName};
%             BinStart=[Info.BinFiles(:).TStart];
%             BinDuration=[Info.BinFiles(:).Duration];
%             HypName={Info.HypnoFiles(:).FileName};
%             HypStart=[Info.HypnoFiles(:).TStart];
%             HypDuration=[Info.HypnoFiles(:).Duration];
% 
%             for nidx=1:length(idx)
% 
%                 Info.BinFiles(nidx).FileName=BinName{idx(nidx)};
%                 Info.BinFiles(nidx).TStart=BinStart(idx(nidx));
%                 Info.BinFiles(nidx).Duration=BinDuration(idx(nidx));
%                 Info.HypnoFiles(nidx).FileName=HypName{idx(nidx)};
%                 Info.HypnoFiles(nidx).TStart=HypStart(idx(nidx));
%                 Info.HypnoFiles(nidx).Duration=HypDuration(idx(nidx));
% 
%             end
%             'exp file loaded'
%         catch ME1
% 
%             ME1.identifier
%             ME1.message
%             for i=1:length(ME1.stack)
%                 ME1.stack(i)
%             end
%             warndlg('exp file not loaded');
%             return;
%         end
%      end
% end
% 
% 


