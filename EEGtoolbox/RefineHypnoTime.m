function RefineHypnoTime(SR1,SR2)
%this function  reshape the hypno if the data hasbeen scored with an error
%in the sampling rate


Info=loadEXP; %load the exp file to process
Duration1=[];
Duration2=[];
for nfile=1:length(Info.BinFiles)
    CurrBinFile=fullfile(Info.FilesDir,Info.BinFiles(nfile).FileName);
    Duration1(nfile)=GetBinDuration(CurrBinFile,Info.NbRecChan)/SR1;
    Duration2(nfile)=GetBinDuration(CurrBinFile,Info.NbRecChan)/SR2;
    
end
figure;
for nhypno=1:length(Info.HypnoFiles)
    CurrHypnoName=fullfile(Info.FilesDir,Info.HypnoFiles(nhypno).FileName);
    fidhypno=fopen(CurrHypnoName);
    OldHypno=fread(fidhypno,'uint16');
    fclose(fidhypno);
    
    %create emptynewhypno
    CurrBinDuration=Duration2(nhypno);
    newHypnoData=zeros(round((CurrBinDuration)),1);
   % formatOut = 'yyyy-mm-ddTHH:MM:SS.FFF';
    %evalute the offset of the hypno starting time
    ts=round(Info.BinFiles(nhypno).TStart*24*3600*1000)/1000; %absolute time in seconde
    OffsetH=5-rem(ts,5); 
    if OffsetH==5
        OffsetH=0;
    end
    
    %resample the hypno;
    newHypnoData=newHypnoData(1:5:end)';
    %remove the first epoch if the offset the more than 2.5
    if OffsetH>2.5
        newHypnoData(1)=[];                
    end

    %check if hypno is longer than the origin file (could be
    %caused by the offset at start
    while ceil(length(newHypnoData)*5+OffsetH)>Duration2(nhypno)
        newHypnoData(end)=[];
    end

    %resample at 1 hz the hypno
    newHypnoData=repmat(newHypnoData,5,1);
    newHypnoData=newHypnoData(:);

    
    
    %fill new hypno with old value at the correct time
  
    for nt=1:length(newHypnoData)
        try
            newHypnoData(nt)=OldHypno(round(nt*SR2/SR1));     
        end
    end
    
    subplot(length(Info.HypnoFiles),1,nhypno);plot(OldHypno,'b');hold on;plot(newHypnoData,'r');
    
    fidhypno=fopen(CurrHypnoName,'w+');
    fwrite(fidhypno,newHypnoData,'uint16');
    fclose(fidhypno);

end

