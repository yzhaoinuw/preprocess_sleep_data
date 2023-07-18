function [FullHypno,TimeScaleAbs,TimeScaleBin,TimeScaleHypno]=ExtractFullHypno(params,CurrExp)
%this function will extract the hypno gram from multiple bin file inluded in the same exp file
%FullHypno  containt the hypnogram
% TimeScaleAbs is the absolute time
% TimeScaleBin is the time in second from the first sample recorded
% TimeScalehypno is the time in second from the first epoch scored


%params containt the file FileInfo whcih is genereted by loadexp  ;
%CurrExp is the experience file number default 1

%exemple to plot the whole hypno

%     params.FileInfo=loadEXP([],'No');
%     [FullHypno,TimeScaleAbs,TimeScaleBin,TimeScaleHypno]=ExtractFullHypno(params,1);
%      figure;plot(TimeScaleHypno,FullHypno


TStart=cat(1,params.FileInfo(CurrExp).HypnoFiles(:).TStart);%*24*3600;
fclose all;
FullHypno=[];

for nHypno=1:length(params.FileInfo(CurrExp).BinFiles);
    fidhyp=fopen(fullfile(params.FileInfo(CurrExp).HypnoFiles(nHypno).Dir,params.FileInfo(CurrExp).HypnoFiles(nHypno).FileName),'r');
    CurrHyp=fread(fidhyp,'uint16');
    idxdeb=etime(datevec(TStart(nHypno)),datevec(TStart(1)))+1;
    idxend=etime(datevec(TStart(nHypno)),datevec(TStart(1)))+length(CurrHyp);
    FullHypno(idxdeb:idxend)=CurrHyp;
    clear('CurrHyp');

end
TimeScaleAbs=linspace(TStart(1),TStart(1)+(length(FullHypno)-1)/3600/24,length(FullHypno));
TimeScaleHypno=0:1:length(FullHypno)-1;

FirstBinTStart=params.FileInfo(CurrExp).BinFiles(1).TStart;
FirstHypnoTStart=params.FileInfo(CurrExp).HypnoFiles(1).TStart;
HypnoOffset=etime(datevec(FirstHypnoTStart),datevec(FirstBinTStart));

TimeScaleBin=TimeScaleHypno+HypnoOffset;