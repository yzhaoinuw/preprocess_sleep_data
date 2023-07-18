function [IdActi]=GetActimetryPosition(CurrExp,CurrBin,handles,params,TimeRel,CurrActi)

%get start time of all video file and the current abs time
TstartCurrBin=params.FileInfo(CurrExp).BinFiles(CurrBin).TStart*24*3600;
Tcurr=TstartCurrBin+TimeRel(1);

allTStartActi=[params.FileInfo(CurrExp).ActiFiles(CurrActi).Files(:).TStart]*24*3600;
allDurationActi=[params.FileInfo(CurrExp).ActiFiles(CurrActi).Files(:).Duration];

% IdActi=0;
% for n=1:size(allTStartActi,2)
%     if Tcurr>=allTStartActi(n) & Tcurr<=allTStartActi+allDurationActi
%     
% end

IdActi=find(Tcurr+60>=allTStartActi & Tcurr<=allTStartActi+allDurationActi,1);

%Sps=params.FileInfo(CurrExp).ActiFiles(CurrActi).Files.SamplingRate;
%Sample=[floor(TimeRel(1).*Sps) ceil(TimeRel(1).*Sps)];