function ExtractChan(Chans,prefix, SubSamp, int16Convert,uint16Convert,params)
%Presicion could be 'uint16' or 'int16'


SampleRate=params.FileInfo.Fs;
Precision=params.FileInfo.Precision;

NewSamplingRate=SampleRate/SubSamp

ChanName=params.FileInfo.ChLabel(Chans)

barre=waitbar(0,'Please wait!');



%create  a low pass filter with a cut frequency half the new sampling rate
 if SubSamp>1
     Filt.hp  = fdesign.lowpass('N,Fc',2,NewSamplingRate/2,SampleRate);
     Filt.Hd = design(Filt.hp,'butter');
 end

%create a new directory to save the new bin files
NewDir=fullfile(params.FilesDir,prefix);
mkdir(NewDir);



 
for nfile=1:length(params.FileInfo.BinFiles)
    
    Filename=fullfile(params.FileInfo.BinFiles(nfile).Dir,params.FileInfo.BinFiles(nfile).FileName);
    
    %process date per hour
    T0Ref =0;
    TmaxRef=GetBinDuration(Filename,params.FileInfo.NbRecChan)/SampleRate;
    TimeLim=T0Ref:3600:TmaxRef;

    if TimeLim~=TmaxRef
        TimeLim=[TimeLim TmaxRef];
    end
    
    %remove short period (less than 60s) and merge the period with the
    %other
    
    idok=diff(TimeLim)>=60;
    TimeLim=TimeLim([true idok]);
    TimeLim(end)=inf;
    
    
    
    Data=[];
    DataOk=[];
    
    Bin2Write=[];
    
    
     %write data in file;

     %keep the initial prefix befor the first _
     pos=strfind(params.FileInfo.BinFiles(nfile).FileName,'_');
     if isempty(pos)==0 && pos(1)>1
         fullprefix=[prefix '_' params.FileInfo.BinFiles(nfile).FileName(1:pos(1)-1)];
     else
         fullprefix=prefix; 
     end
    currfilename=fullfile(NewDir,[fullprefix '_'  params.FileInfo.BinFiles(nfile).FileName(end-26:end)]);
    fidnew=fopen(currfilename,'w+');
    
    for nhour=1:length(TimeLim)-1
        
        waitbar(nhour/(length(TimeLim)-1)/length(params.FileInfo.BinFiles)+(nfile-1)/length(params.FileInfo.BinFiles),barre);
        
        
         T0=TimeLim(nhour);
         Tmax=TimeLim(nhour+1);
        
        

            fid(nfile)=fopen(Filename,'r+');
            %[Data]=ReadBinFile(SampleRate,0,inf,fid(nfile),params.FileInfo.NbRecChan,1,Precision);
            [Data]=ReadBinFile(SampleRate,T0,Tmax,fid(nfile),params.FileInfo.NbRecChan,1,Precision);

            DataOk=Data(Chans,:);

            clear Data;

            fclose(fid(nfile));


            %FilterData
            NewFileData=[];
            Bin2Write;

            if SubSamp>1
                %add a low pass filter at half the new sampling rate
                NewFileData=[filtfilt(Filt.Hd.sosMatrix,Filt.Hd.ScaleValues,DataOk')]';
                NewFileData=NewFileData(:,1:SubSamp:end);
                %NewFileData=[resample(DataOk',NewSamplingRate,SampleRate)]';

            else
                NewFileData=[DataOk];
            end

            NewPrecision= Precision;
            if uint16Convert==1 && strcmpi(Precision,'uint16')==0
                NewFileData=NewFileData+2^15;
                NewPrecision='uint16';

            end

            if int16Convert==1 && strcmpi(Precision,'int16')==0
                NewFileData=NewFileData-2^15;
                 NewPrecision='int16';
            end
            
            %Bin2Write=eval(sprintf('cat(2,Bin2Write,%s(NewFileData));',NewPrecision));
            
            Bin2Write=eval(sprintf('%s(NewFileData);',NewPrecision));
         
            if uint16Convert==1
                fwrite(fidnew,Bin2Write,'uint16');
            elseif int16Convert==1 
                 fwrite(fidnew,Bin2Write,'int16');
            end
        
    end

   

    
    fclose(fidnew);

end

%create a file information about the conversion
fidInfo=fopen(fullfile(NewDir,[prefix 'FileInfo.txt']),'W+');

fprintf(fidInfo,'Conversion information\n\n');
fprintf(fidInfo,'Sampling Rate: \t %f \n',NewSamplingRate);
fprintf(fidInfo,'Nb Chan:\t %f \n',length(ChanName));
fprintf(fidInfo,'Precision:\t %s \n',NewPrecision);

fprintf(fidInfo,'Gain:\t [');
for n=1:length(ChanName)
    fprintf(fidInfo,' %f',params.FileInfo.Gain(Chans(n)));
end
fprintf(fidInfo,' ]\n');

fprintf(fidInfo,'MaxRange:\t [');
for n=1:length(ChanName)
    fprintf(fidInfo,' %f',params.FileInfo.Range(Chans(n)));
end
fprintf(fidInfo,']\n');

NewOffset=params.FileInfo.Offset(Chans);

if uint16Convert==1 && strcmp(Precision,'uint16')==0
    %conversion from int16 to uint16
    NewOffset=params.FileInfo.Offset(Chans)+params.FileInfo.Range(Chans)./2;

elseif int16Convert==1 && strcmp(Precision,'int16')==0
    %conversion from uint16 to int16
    NewOffset=params.FileInfo.Offset(Chans)-params.FileInfo.Range(Chans)./2;

else
    NewOffset=params.FileInfo.Offset(Chans);

end


fprintf(fidInfo,'Offset:\t [');
for n=1:length(NewOffset)
    fprintf(fidInfo,' %f',NewOffset(n));
end
fprintf(fidInfo,']\n');


fprintf(fidInfo,'ChanName:\t {');
for n=1:length(ChanName)
    fprintf(fidInfo,' ''%s''',ChanName{n});
end
fprintf(fidInfo,'}\n');
fclose(fidInfo);








close(barre);




warndlg('done');
    
