function Quantif=StateQuantification(Info,T0,Duration,StateCycle,MinEpisodDuration,SeqBinSize,State2Process,MainOutDir,SaveOn)
%This script process the quantification analysis from an exp file
%info the the struct generatated by loadEXP
%T0  TStart from the fisrt bin file in the exp in min
%Duration to analyse in min
%StateCycle is a matrix that contain the code which are considerg as a sleep episod
%MinEpisodDuration minimal duration to consider in s
%SeqBinSize in min duration of the sequential step
%State2Process state to process for quantification
%SaveOn =1 if data saved


%get the exp if not present

if isempty(Info)
    try
        load tempdir
    catch
        DefaultDir=cd;
    end
    [f,d]=uigetfile('*.exp','pick the exp file to process',DefaultDir);
    
    CurrExpFile=fullfile(d,f);
    Info=loadEXP(CurrExpFile,'no');
else
    d=Info.FilesDir;
    f=Info.Filename;
end
params.FileInfo=Info;
Header.Info=Info;

[FullHypno,TimeScaleAbs,TimeScaleBin,TimeScaleHypno]=ExtractFullHypno(params,1);


if T0+Duration>TimeScaleHypno(end)/60
   Duration = TimeScaleHypno(end)/60-T0; 
end


% %Extract time lim from hypno
[Header,Data]=ExtractTimeLimFromHypno(Header,MinEpisodDuration,T0,Duration,[],[],[],[],0);


Quantif.Settings.T0=T0;
Quantif.Settings.Duration=Duration;
Quantif.Settings.SeqBinSize=SeqBinSize;
Quantif.Settings.State2Process=State2Process;
idssArtefNd=~(Data.Event.Hypno==0 | Data.Event.Hypno==5); %State 5 (Undetermined) & 0 (Artefact) are 
                                                          %substract to the total duration for percentage calcultation
                                                          
SeqLim=T0:SeqBinSize:T0+Duration; 

if   SeqLim(end)~=T0+Duration
	SeqLim=[SeqLim T0+Duration];
end
                                                          
                                                          
for n=1:length(State2Process)
    


    %trouver les index des episodes a traiter
    id=Data.Event.Hypno==State2Process(n);


    Quantif.State(n).TotDuration=sum(Data.Event.Duration(id))/60;
    Quantif.State(n).PerDuration=Quantif.State(n).TotDuration/(sum(Data.Event.Duration(idssArtefNd))/60)*100;
    Quantif.State(n).nbEpisod=sum(id);
    Quantif.State(n).MeanEpisodDuration=mean(Data.Event.Duration(id))/60;
    Quantif.State(n).StdEpisodDuration=std(Data.Event.Duration(id))/60;
    Quantif.State(n).MedianEpisodDuration=median(Data.Event.Duration(id))/60;
    StatePos=[Header.Info.State.Code]==State2Process(n);
    Quantif.State(n).Label=Header.Info.State(StatePos).Label;
    Quantif.State(n).Code=Header.Info.State(StatePos).Code;
    Quantif.State(n).Color=Header.Info.State(StatePos).Color;
     
    for Nseq=1:length(SeqLim)-1
        IdSeqok=id & Data.Event.Deb/60>=SeqLim(Nseq) &  Data.Event.Deb/60<=SeqLim(Nseq+1);
        
        Quantif.State(n).SeqTotDuration(Nseq)=sum(Data.Event.Duration(IdSeqok))/60;
        Quantif.State(n).SeqPerDuration(Nseq)=Quantif.State(n).TotDuration/(sum(Data.Event.Duration(idssArtefNd))/60)*100;
        Quantif.State(n).SeqnbEpisod(Nseq)=sum(IdSeqok);
        Quantif.State(n).SeqMeanEpisodDuration(Nseq)=mean(Data.Event.Duration(IdSeqok))/60;
        Quantif.State(n).SeqStdEpisodDuration(Nseq)=std(Data.Event.Duration(IdSeqok))/60;
        Quantif.State(n).SeqMedianEpisodDuration(Nseq)=median(Data.Event.Duration(IdSeqok))/60;
        Quantif.State(n).SeqLimDeb(Nseq)=SeqLim(Nseq);
        Quantif.State(n).SeqLimEnd(Nseq)=SeqLim(Nseq+1);
    end
    
    
end










if isempty(StateCycle)==0
    %extract cycle
    %sort the data from the state code deine by StateCycle
    %Extract time lim from hypno
    [~,DataCycle]=ExtractTimeLimFromHypno(Header,MinEpisodDuration,T0,Duration,[],[],StateCycle,[],0);

    %find transition state 1 vers 2

    idS1=DataCycle.Event.Hypno==StateCycle(1);
    idS2=DataCycle.Event.Hypno==StateCycle(2);

    H=DataCycle.Event.Hypno;
    H(idS1)=1;
    H(idS2)=2;


    idOk=[DataCycle.Event.End(1:end-1)==DataCycle.Event.Deb(2:end) 0] & idS1;

    %merge cycle
    DataCycleMerge.Deb=DataCycle.Event.Deb(find(idOk));
    DataCycleMerge.End=DataCycle.Event.End(find(idOk)+1);
    DataCycleMerge.Duration=DataCycleMerge.End-DataCycleMerge.Deb;
    DataCycleMerge.DurationState1=DataCycle.Event.Duration(find(idOk));
    DataCycleMerge.DurationState2=DataCycle.Event.Duration(find(idOk)+1);

    %outname

    
    Quantif.DataCycleMerge=DataCycleMerge;
else
    Quantif.DataCycleMerge=[];
end

if SaveOn==1

    %save in a mat file with the same name than the exp file with the suffix 'Quantif'
    if exist('MainOutDir')==0 || (exist('MainOutDir')~=0 && isempty(MainOutDir)==1)
        CurrDir=Info.FilesDir;
    else
        CurrDir=MainOutDir;
    end

    mkdir(sprintf('%s',CurrDir));
    MatFilename=fullfile(CurrDir,['Quantif_' Header.Info.ExpFileName(1:end-4) '.mat']);
    save(MatFilename,'Quantif');
    
    
    %save in xls
    XlsxFileName=fullfile(CurrDir,['Quantif_' Header.Info.ExpFileName(1:end-4) '.xlsx']);
    xlswrite(XlsxFileName,{'State' 'Duration (min)' 'Duration (%)' 'nb episod' 'mean episod duration (min)' 'StdEpisodDuration (min)' 'MedianEpisodDuration (min)'},'Quantif');

   T={};
   T{1}='State';
   
   for nlim=1:length(Quantif.State(1).SeqLimDeb)
      T{nlim+1}=sprintf('%5.1f',Quantif.State(1).SeqLimDeb(nlim));
   end
   T{nlim+2}=sprintf('%5.1f',Quantif.State(1).SeqLimEnd(nlim));
   xlswrite(XlsxFileName,T,'Sequential TotDuration');
   
  Tab=[];
  TabSeqTotDuration=[];
   for nstate=1:length(State2Process)
       Tab(nstate,1:7)=[Quantif.State(nstate).Code Quantif.State(nstate).TotDuration Quantif.State(nstate).PerDuration ...
          Quantif.State(nstate).nbEpisod  Quantif.State(nstate).MeanEpisodDuration Quantif.State(nstate).StdEpisodDuration Quantif.State(nstate).MedianEpisodDuration];
      
       TabSeqTotDuration(nstate,1:length(Quantif.State(nstate).SeqLimDeb)+1)=[Quantif.State(nstate).Code Quantif.State(nstate).SeqTotDuration];
       
   
   end
   xlswrite(XlsxFileName,Tab,'Quantif','A2');
   xlswrite(XlsxFileName,TabSeqTotDuration,'Sequential TotDuration','A2');
    
    
    
    if isempty(StateCycle)==0
       OutNamePrefix=['Quant-' f(1:end-4)];

        C=cell(length(DataCycleMerge.Deb)+1,6);
        Data=num2cell([DataCycleMerge.Deb' DataCycleMerge.End' DataCycleMerge.Duration']);
        C(2:end,1:3)=Data;
        C(1,1:6)={'Start' 'End' 'Duration (s)' 'MeanDuration (s)' 'Std Duration(s)' 'Number'};
        stat=[mean(DataCycleMerge.Duration) std(DataCycleMerge.Duration) length(DataCycleMerge.Duration)];
        C(2,4:end)=num2cell(stat);

        xlswrite(XlsxFileName,C,'StateCycle');
        save(MatFilename,'DataCycleMerge','-append');
        
 
    end
    
    
    

end

