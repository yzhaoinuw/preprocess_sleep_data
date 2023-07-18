function ChangeState(Info,T0,TEnd,StateCode1,StateCode2,NewStateCode,StateMaxDuration,NewHypnoPref)

%T0,TEnd are the start time and end time in min
%this function rescore State1 in NewState if State1 is surrounded by
%State2 and if the episod last less than StateMaxDuration (in s), State2 could be empty in that case all state lasting less than
%StateMaxDuration will be changed
%NewHypnoPref is the the exp file



%% open exp file
if isempty(Info)==1
    Info=loadEXP([],'No');
end
Directory=Info.FilesDir;


%get original hypno
params.FileInfo=Info;
[FullHypno,TimeScaleAbs,TimeScaleBin,TimeScaleHypno]=ExtractFullHypno(params,1);

nbEpisodChanged=0;
if isempty(NewHypnoPref)==0 %erase the original hypno
    %Duplicate the hypno and EXP;
    OldExpFilenameExt=Info.ExpFileName;
    NewEXPFilename=fullfile(Directory,[NewHypnoPref '.exp']);

    [PATHSTR,NewEXPFilename,EXT] = fileparts(NewEXPFilename);
    NewFile=Info;
    NewEXPFilenameExp=[NewEXPFilename '.exp'];

    formatIn = 'yyyy-mm-dd_HH-MM-SS-FFF';

    BoutsDuration=1;

    %load old exp struct
    [ s ] = Xml2struct(fullfile(Directory,OldExpFilenameExt));

    %update the fileorigin in the exp struct
    s.Animal.FileOrigin.Text=Directory;
    NewFile.Info=Directory;

    PrefixeFolder=['Hypno_' NewEXPFilename];

    mkdir(fullfile(Directory,PrefixeFolder));


    for nHyp=1:length(NewFile.HypnoFiles)
        %update the structure of the nex exp file
        NewFile.HypnoFiles(nHyp).FileName=fullfile(PrefixeFolder,[NewEXPFilename '_' Info.HypnoFiles(nHyp).FileName(end-length(formatIn)-1:end)]); 

        if exist(fullfile(Directory,NewFile.HypnoFiles(nHyp).FileName))==1
            Ans=questdlg('File already exist, would you like to replace it?');
            if strcmp(Ans,'Yes')    
                %duplicate hypnofile in a new file
                copyfile(fullfile(Info.HypnoFiles(nHyp).Dir,Info.HypnoFiles(nHyp).FileName),...
                    fullfile(Directory,NewFile.HypnoFiles(nHyp).FileName));
            end
        else

             copyfile(fullfile(Info.HypnoFiles(nHyp).Dir,Info.HypnoFiles(nHyp).FileName),...
                    fullfile(Directory,NewFile.HypnoFiles(nHyp).FileName));
        end

        %update the new struct
        s.Animal.Hypnogram.Files.File{nHyp}.FileName.Text=NewFile.HypnoFiles(nHyp).FileName;

    end


    %duplicate the matfile with the new name    
    copyfile(fullfile(Directory,[OldExpFilenameExt(1:end-3) 'mat']),...
            fullfile(Directory,[NewEXPFilename '.mat']));


    %save the nex EXP;
    XMLCurrFileName=[NewEXPFilename '.xml'];
    struct2xml(s,fullfile(Directory,XMLCurrFileName));
    copyfile(fullfile(Directory,XMLCurrFileName),fullfile(Directory,NewEXPFilenameExp),'f');
    delete(fullfile(Directory,XMLCurrFileName));
    clear Info NewFile;


    %reload the newexpfile;
    MaxRangeType='No';

    Info=loadEXP(fullfile(Directory,NewEXPFilenameExp),MaxRangeType);
    
end
%For each hypno

NbHypno=length(Info.HypnoFiles);

for nHypno=1:NbHypno


    %open bin and hypno file
    %fidbin=fopen(fullfile(params.FilesDir,Info.BinFiles(nHypno).FileName),'r');
    fidhyp=fopen(fullfile(Info.HypnoFiles(nHypno).Dir,Info.HypnoFiles(nHypno).FileName),'r+');
    
    %get the offset between hypno an bin
    offsetHypBin=etime(datevec(Info.HypnoFiles(nHypno).TStart),datevec(Info.BinFiles(nHypno).TStart));
    CurrHypno=fread(fidhyp,'uint16')';
    
    tStartHypno=offsetHypBin+etime(datevec(Info.BinFiles(nHypno).TStart),datevec(Info.BinFiles(1).TStart));
    tEndHypno=tStartHypno+length(CurrHypno);

  
    %Generate time matrix
    Thypno=linspace(tStartHypno,tEndHypno,length(CurrHypno));
    AbsPos=1:length(CurrHypno);
    %replace data if needed
    if sum(Thypno>=T0 &  Thypno<=TEnd)>1
        
        %Detect all episod of Stats1
        idState1=CurrHypno==StateCode1;
        
        %detect all deb and end of each episod
            % %extract start and end times
        Diff1=diff([0 idState1 0])==1;
        Diff1_1=diff([0 idState1 0])==-1;
        State1.Deb=AbsPos(Diff1(1:end-1));
        State1.End=AbsPos(Diff1_1(2:end));
        
        %remove all Episod longer than StateMaxDuration
        
        idok=(State1.End-State1.Deb)<=StateMaxDuration;
        State1.Deb=State1.Deb(idok);
        State1.End=State1.End(idok);
        State1.Duration=State1.End-State1.Deb;
        
        %get the state befor and after the episod
        CurrHypnoPrev=[0 CurrHypno];%add a zero to avoir any pb
        CurrHypnAfter=[CurrHypno 0];%add a zero to avoir any pb
        State1.PrevState=CurrHypnoPrev(State1.Deb);
        State1.AftState=CurrHypnAfter(State1.End+1);
        
       
        
        
        %for each episod
        for n=1:length(State1.Deb)
            
            AbsTimeDeb=Thypno(State1.Deb(n));
            AbsTimeEnd=Thypno(State1.End(n));
            
            if State1.PrevState(n)==StateCode2 && State1.AftState(n)==StateCode2  && AbsTimeDeb/60>=T0 && AbsTimeEnd/60<=TEnd
                CurrHypno(State1.Deb(n):State1.End(n))=NewStateCode;
                nbEpisodChanged=nbEpisodChanged+1;
            end
            
            
        end


        %write new hypno
        Data2Write=uint16(CurrHypno);
        fseek(fidhyp,0,'bof');
        fwrite(fidhyp,Data2Write,'uint16');
    end
    fclose(fidhyp);

 
end
nbEpisodChanged
'LM scoring done'




