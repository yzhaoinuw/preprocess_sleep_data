function [ModeState,AllCode]=GetStateBetweenLim(Info,T0,Duration)
%this function gives the mode state and all the code between limit
%T0 in s , time relative to start
%Duration in s duration to extract

params.FileInfo=Info;

[FullHypno,~,TimeScaleBin,~]=ExtractFullHypno(params,1);

AllCode=FullHypno(TimeScaleBin>=T0 &  TimeScaleBin<=T0+Duration);
ModeState=mode(AllCode);

