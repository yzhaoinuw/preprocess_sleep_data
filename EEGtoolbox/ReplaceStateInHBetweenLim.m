
function ReplaceStateInHBetweenLim(Info,NewStateCode,T0,TEnd)
%StateCode code to force in H
%T0 in s from first bin file start
%TEnd in s from first bin file start

nbBinFile=length(Info.BinFiles);

for nbin=1:nbBinFile


    %open bin and hypno file
    %fidbin=fopen(fullfile(Directory,Info.BinFiles(nbin).FileName),'r');
    fidhyp=fopen(fullfile(Info.HypnoFiles(nbin).Dir,Info.HypnoFiles(nbin).FileName),'r+');
    
    CurrHypno=fread(fidhyp,'uint16');
    NbSample=length(CurrHypno);
    tStartHypno=etime(datevec(Info.HypnoFiles(nbin).TStart),datevec(Info.BinFiles(1).TStart));
    tEndHypno=tStartHypno+length(CurrHypno)-1;

    if isinf(TEnd)==1
        TEndCurr=tEndHypno;
    else
        TEndCurr=TEnd;
    end

    %Generate time matrix
    Thypno=linspace(tStartHypno,tEndHypno,NbSample);
    
    
    T0SecBefAbs=floor(Info.BinFiles(1).TStart*24*3600+T0)/24/3600;
    T0SecBefRel=etime(datevec(T0SecBefAbs),datevec(Info.BinFiles(1).TStart));
    TEndSecAFterAbs=ceil(Info.BinFiles(1).TStart*24*3600+TEndCurr)/24/3600;
    TEndSecAFter=etime(datevec(TEndSecAFterAbs),datevec(Info.BinFiles(1).TStart));
    
    %replace data if needed
    if sum(Thypno>=T0SecBefRel &  Thypno<=TEndSecAFter)>=1
        id=Thypno>=T0SecBefRel &  Thypno<TEndSecAFter;

        if length(NewStateCode)==1
            CurrHypno(id)=NewStateCode;
        else
            
            TnewCode=T0:length(NewStateCode)+T0;
            
            first=find(TnewCode>=Thypno(1),1,'first');
            last=find(TnewCode<=Thypno(end),1,'last');
            FirstCurrHypPos=find(Thypno>=TnewCode(first),1,'first');
            shortNewH=NewStateCode(first:last-1);
          
            
            CurrHypno(FirstCurrHypPos:FirstCurrHypPos+length(shortNewH)-1)=shortNewH;
           % id=Thypno>=T0 &&  TnewCode<=tEndHypno;
            
            
            
        end

        %write new hypno
        Data2Write=uint16(CurrHypno);
        fseek(fidhyp,0,'bof');
        fwrite(fidhyp,Data2Write,'uint16');
    end
    fclose(fidhyp);

 
end