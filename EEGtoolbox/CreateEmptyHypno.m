function [DataH,HStart]=CreateEmptyHypno(BinDuration, BinStart,EpochDuration)
%this function create an empty hypno starting from the first multiple of
%the epoch duration

%BinDuration in s, the duration of the binary file
%BinStart, time start of the bin file in vecnum
%EpochDuration in s is the epoch duration

%DataH is a matrix containing the hypnocode at Hz in uint16
%HStart, time start of the H file in datenum


%evalute the offset of the hypno starting time
BinStartVec=datevec(BinStart); %time in vec
% formatOut = 'yyyy-mm-dd HH:MM:SS.FFF'
% datestr(BinStart,formatOut)
BinStartSec=BinStartVec(6);
OffsetStart=EpochDuration-rem(BinStartSec,EpochDuration);
if OffsetStart==EpochDuration
    OffsetStart=0;
end

HStart=(BinStart*3600*24+OffsetStart)/3600/24;


% datestr(HStart,formatOut)

BinEnd=(BinStart*24*3600+BinDuration)/24/3600;
BinEndVec=datevec(BinEnd);
BinEndSec=BinEndVec(6);
%datestr(BinEnd,formatOut)
OffsetEnd=rem(BinEndSec,EpochDuration);

HEnd=(BinStart*3600*24-OffsetEnd+BinDuration)/3600/24; 
% datestr(HEnd,formatOut)
nbsamp=floor(etime(datevec(HEnd),datevec(HStart)));

DataH=zeros(nbsamp,1);
