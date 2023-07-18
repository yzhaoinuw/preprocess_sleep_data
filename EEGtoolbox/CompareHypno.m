function [DataComp]=CompareHypno(StateToReplace,State2Evaluate,T0,TEnd,FileNameRef,FileNameScorageAuto,OutDir)
%StateToReplace is a matrix with 2 line and n row, each raw containt the orginal state and the new
%State2Evaluate comtain the state code to analyse;
%T0 and TEnd are in Min

%exemple [CP,DataComp]=CompareHypno([1 2 3 4 7 8 9;1 1 3 1 7 7 1],[1 3 7],0,1440)

if nargin<=4
    FileNameRef=[];
    FileNameScorageAuto=[];
elseif nargin<=5
    FileNameScorageAuto=[];
end


%load hypno 1 ref hypno
' reference hypno '
if isempty(FileNameRef)==1
    Info1=loadEXP([],'No');
    params1.FileInfo=Info1;
    [Hypno1,TimeScaleAbs1,TimeScaleBin1,TimeScaleHypno1]=ExtractFullHypno(params1,1);
else
    Info1=loadEXP(FileNameRef,'No');
    params1.FileInfo=Info1;
    [Hypno1,TimeScaleAbs1,TimeScaleBin1,TimeScaleHypno1]=ExtractFullHypno(params1,1);
end

if nargin<=6
    OutDir=Info1.FilesDir;
end

%load Hypno2
if isempty(FileNameScorageAuto)==1
    Info2=loadEXP([],'No');
    params2.FileInfo=Info2;
    [Hypno2,TimeScaleAbs2,TimeScaleBin2,TimeScaleHypno2]=ExtractFullHypno(params2,1);
else
    Info2=loadEXP(FileNameScorageAuto,'No');
    params2.FileInfo=Info2;
    [Hypno2,TimeScaleAbs2,TimeScaleBin2,TimeScaleHypno2]=ExtractFullHypno(params2,1);
end

%Check if AbsTimeStart are the same

if TimeScaleAbs1(1)==TimeScaleAbs2(1)

    %reduce Hypnoduration
    if nargin<=2
        T0=0;
        TEnd=inf;
    end

    
    
    
    id1=TimeScaleHypno1>=T0*60 & TimeScaleHypno1<=TEnd*60;
    id2=TimeScaleHypno2>=T0*60 & TimeScaleHypno2<=TEnd*60;
    
    Hypno1=Hypno1(id1);
    Hypno2=Hypno2(id2);


    %replace code
    for n=1:size(StateToReplace,2)
        id=Hypno1==StateToReplace(1,n);
        Hypno1(id)=StateToReplace(2,n);
        id=Hypno2==StateToReplace(1,n);
        Hypno2(id)=StateToReplace(2,n);
    end
    
        
    if exist('State2Evaluate')
        id1H=ismember(Hypno1,State2Evaluate);
        id2H=ismember(Hypno2,State2Evaluate);
        
        Hypno1=Hypno1(id1H & id2H);
        Hypno2=Hypno2(id1H & id2H);
        TimeScaleHypno1=TimeScaleHypno1(id1H & id2H);
        TimeScaleHypno2=TimeScaleHypno2(id1H & id2H);
        
    end
    
    
    
     %plot the 2 hypno
    fcomp=figure;
    h(1)=subplot(2,1,1);
    plot(TimeScaleHypno1/60,Hypno1);
    ylabel('HypnoRef');
    grid on;
    h(2)=subplot(2,1,2);
    plot(TimeScaleHypno2/60,Hypno2);
    linkaxes(h,'xy')
    xlabel('Time (min)')
    ylabel('Hypno2');
    grid on;
    title(Info2.Name);

    CompFigFilename=fullfile(OutDir,sprintf('Comp %s vs %s .fig',Info1.ExpFileName(1:end-4),Info2.ExpFileName(1:end-4)));
    CompbmpFilename=fullfile(OutDir,sprintf('Comp %s vs %s .bmp',Info1.ExpFileName(1:end-4),Info2.ExpFileName(1:end-4)));
    saveas(fcomp,CompFigFilename);
    print(fcomp,CompbmpFilename,'-dbmp','-r600');



    %compare hypno
%      CP = classperf(Hypno1,Hypno2);
%      get(CP)
     DataComp.AllCode1=[Info1.State.Code];
     DataComp.AllName1={Info1.State.Label};
     DataComp.AllCode2=[Info2.State.Code];
     DataComp.AllName2={Info2.State.Label};
     
     StateCode=unique(Hypno1);
     
     DataComp.StateName1=DataComp.AllName1(ismember(DataComp.AllCode1,StateCode));
     DataComp.StateName2=DataComp.AllName2(ismember(DataComp.AllCode2,StateCode));
     %StateName=

     DataComp.ConfMatrix=cell(length(DataComp.StateName1)+1,length(DataComp.StateName1)+1);
     DataComp.ConfMatrix(2:end,1)=DataComp.StateName1;
     DataComp.ConfMatrix(1,2:end)=DataComp.StateName1';
     %C=get(CP,'CountingMatrix');
     
    % DataComp.ConfMatrix(2:end,2:end)=num2cell(C(1:end-1,:));
     C=confusionmat(Hypno1,Hypno2)';
     DataComp.ConfMatrix(2:end,2:end)=num2cell(C);
     ErrorDistributionByClass=sum(C.*~eye(size(C,1)),1)';
     SampleDistributionByClass=sum(C,2);
     errorrate=sum(ErrorDistributionByClass)/sum(SampleDistributionByClass);
     CorrectRate=1-errorrate;

     DataComp.CorrectRate=CorrectRate;
     
     DataComp.CorrectRatePerState=cat(2,[DataComp.StateName1';'Total'],[num2cell(SampleDistributionByClass./(SampleDistributionByClass+ErrorDistributionByClass));{CorrectRate}]);
     for nstate=1:size(C,1)
         TP(nstate)=C(nstate,nstate);
         FP(nstate)=sum(C(nstate,:))-TP(nstate);
         Idnot=true(size(C,1),size(C,2));
         Idnot(nstate,:)=false;
         Idnot(:,nstate)=false;
         TN(nstate)=sum(C(Idnot));
         FN(nstate)=sum(C(:,nstate))-TP(nstate);
     end
     
     DataComp.Sensitivity=TP./(TP+FN);
     DataComp.Specificity=TN./(FP+TN);
     DataComp.PositivePredictiveValue=TP./(TP+FP);
     DataComp.NegativePredictiveValue=TN./(FN+TN);
     
     DataComp.MeanSpecificity=mean(DataComp.Specificity);
     DataComp.MeanSensitivity=mean(DataComp.Sensitivity);
     DataComp.MeanPositivePredictiveValue=mean(DataComp.PositivePredictiveValue);
     DataComp.MeanNegativePredictiveValue=mean(DataComp.NegativePredictiveValue);
     
%      DataComp.CorrectRate
%      DataComp.Specificity=get(CP,'Specificity');
%      DataComp.Sensitivity=get(CP,'Sensitivity');
%      DataComp.CorrectRate=get(CP,'CorrectRate');
%    DataComp.CorrectRatePerState=cat(2,[DataComp.StateName1';'Total'],[num2cell(get(CP,'SampleDistributionByClass')./(get(CP,'SampleDistributionByClass')+get(CP,'ErrorDistributionByClass')));{DataComp.CorrectRate}]);
%      DataComp.PositivePredictiveValue=get(CP,'PositivePredictiveValue');
%      DataComp.NegativePredictiveValue=get(CP,'NegativePredictiveValue');
try
     DataComp.kappa = cohenskappa(Hypno1,Hypno2);
end
%      DataComp.PositiveLikelihood=get(CP,'PositiveLikelihood');
%      DataComp.CP=CP;
     DataComp
     
     DataMatName=fullfile(OutDir,sprintf('Comp %s vs %s .mat',Info1.ExpFileName(1:end-4),Info2.ExpFileName(1:end-4)));
     save(DataMatName,'DataComp');
else
    'Abs Time Start are not the same'
end