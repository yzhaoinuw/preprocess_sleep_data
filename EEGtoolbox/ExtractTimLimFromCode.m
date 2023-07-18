function [TimDeb,TimeEnd]=ExtractTimLimFromCode(Data,RawTime,Timewin,Code,MinDuration,MinSpaceBetween2Lim,MinDuration2)
%Data is a line vector containing the data to threshold
%Raw time is the time vector of Data in s, Raw time is the middle of each windows used for computing Data
%Timewin is the size of the windows used for computing each value in Data(in s)
%Code 2 Code of the hypnogram to select
%MinDuration is the min event duration in s all events lasting less to MinDuration are removed before any other operation
%MinSpaceBetween2Lim is the minimum intervalle without merging time limit All gap shorter to this value is filled
%MinDuration2 is the min event duration in s all events lasting less to MinDuration2 are removed done after Filling Gap


%Get Code index
Id1=Data==Code;
 
% %extract start and end times
Diff1=diff([0 Id1 0])==1;
Diff1_1=diff([0 Id1 0])==-1;
TimDeb=RawTime(Diff1(1:end-1))-Timewin/2;
TimeEnd=RawTime(Diff1_1(2:end))+Timewin/2;
Duration1=TimeEnd-TimDeb;

%remove shorter Oscillation phase1
Id2Short1=Duration1<MinDuration;
TimDeb(Id2Short1)=[];
TimeEnd(Id2Short1)=[];

%merge event with to short interval between them 
IOI1=TimDeb(2:end)-TimeEnd(1:end-1);

ShortIOIIdStart=[true IOI1>MinSpaceBetween2Lim];
ShortIOIIdEnd=[IOI1>MinSpaceBetween2Lim true];

TimDeb=TimDeb(ShortIOIIdStart);
TimeEnd=TimeEnd(ShortIOIIdEnd);
Duration1=TimeEnd-TimDeb;

%remove shorter Oscillation phase2
Id2Short1=Duration1<MinDuration2;
TimDeb(Id2Short1)=[];
TimeEnd(Id2Short1)=[];