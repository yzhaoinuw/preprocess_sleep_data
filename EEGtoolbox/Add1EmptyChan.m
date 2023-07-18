function Add1EmptyChan(Pos,BinVal,Prefix)
%add 1 channel at the Pos position with the BinVal add the prefix Prefix to the file name

Info=loadEXP;
CurrDir=Info.FilesDir;
nbchan=Info.NbRecChan;
Precision=Info.Precision;

for nbbinfile=1:length(Info.BinFiles)
    Filename=Info.BinFiles(nbbinfile).FileName;
    fid=fopen([CurrDir Filename],'r+');
    Data=[];
    Data=fread(fid,[nbchan inf],Precision);
    fclose(fid);
    
    Col=repmat(BinVal,1,size(Data,2));
    if Pos==1
        Data=[Col;Data];
    elseif Pos==nbchan+1
        Data=[Data;Col];
    else
        Data=[Data(1:Pos-1,:);Col;Data(Pos+1:end,:)];
    end
    
    
     fidnew=fopen([CurrDir Prefix Filename],'w+');
     fwrite(fidnew,Data,Precision);
     fclose(fidnew);
    
    
end