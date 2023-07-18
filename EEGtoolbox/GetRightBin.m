function [GoodBin,ClosestBin]=GetRightBin(Info,Time,TimeMode)

    %Timemode is 1 if absolute Time format datenum
    %timemode is 2 if relative mode (time from start in secondes)
    %ClosestBin is the previous and next bin file found if the Time is
    %between 2 bin files

    %generate the timebase per second to search the correct bin    
    GoodBin=NaN;
 
    DurationBinFile=cat(1,Info.BinFiles(:).Duration);
    


%     TStart=cat(1,Info.BinFiles(:).TStart)*24*3600;

%     if TimeMode==1 
%         for nFile=1:length(Info.BinFiles)
%             TimeFileDebAbs=TStart(nFile);
%             TimeFileEndAbs=TStart(nFile)+DurationBinFile(nFile);
%             Time=Time*24*3600;
%             if Time>=TimeFileDebAbs && Time<TimeFileEndAbs
%                 GoodBin=nFile;
%             end
%         end
%     end
% 
% 
%     if TimeMode==2
%         for nFile=1:length(Info.BinFiles)
%             TimeFileDebRel=TStart(nFile)-TStart(1);
%             TimeFileEndRel=TStart(nFile)-TStart(1)+DurationBinFile(nFile);
% 
%             if Time>=TimeFileDebRel && Time<TimeFileEndRel
%                 GoodBin=nFile;
%             end
%         end
%     end
TStart=cat(1,Info.BinFiles(:).TStart);
if TimeMode==1 
    CurrAbsTime=Time;
elseif TimeMode==2 
	CurrAbsTime=TStart(1)+Time/24/3600;
end


T2DebFiles=etime(datevec(repmat(CurrAbsTime,size(TStart))),datevec(TStart(:)));
T2EndFiles=etime(datevec(TStart(:)+DurationBinFile(:)/24/3600),datevec(repmat(CurrAbsTime,size(TStart))));

GoodBin=find(T2DebFiles>=0 & T2EndFiles>=0);
ClosestBin=[GoodBin GoodBin];

if isempty(GoodBin)==1
    GoodBin=NaN;
    ClosestBin=[find(T2DebFiles>=0,1,'last') find(T2EndFiles>=0,1,'first')];
    'no bin file found for this limit'
elseif length(GoodBin)>=2%if bin superpositon
    GoodBin=GoodBin(1);%take the first bin that include the time requested
end

end
