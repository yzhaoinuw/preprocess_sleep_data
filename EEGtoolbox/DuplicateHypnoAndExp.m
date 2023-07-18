function [NewFileInfo,Status]=DuplicateHypnoAndExp(paramsOldEXP,NewEXPFilename)
%this function save the current exp file and hypno in a new folder by
%generating hypno with prefix enad exp with prefix
%paramsOldEXP containt the params with the fileinfo generated from loadEXP
%in a structure like params.FileInfo=loadEXP
%NewEXPFilename is the new name of the hypno and exp file without extention
%NewFileInfo containt  the exp information of the new exp file

 try
    OldExpFilenameExt=paramsOldEXP.FileInfo.Filename;
    NewEXPFilenameExp=[NewEXPFilename '.exp'];
    formatIn = 'yyyy-mm-dd_HH-MM-SS-FFF';
    %load old exp struct
    [ s ] = Xml2struct(fullfile(paramsOldEXP.FileInfo.FilesDir,OldExpFilenameExt));

    %update the fileorigin in the exp struct
    s.Animal.FileOrigin.Text=paramsOldEXP.FileInfo.FilesDir;

    NewFileInfo=paramsOldEXP.FileInfo;
    NewFileInfo.Filename=NewEXPFilenameExp;
    NewFileInfo.ExpFileName=NewEXPFilenameExp;
    NewFileInfo.FilesDir=paramsOldEXP.FileInfo.FilesDir;
%     NewFileInfo.FilesDir
    %NewFileInfo.FileInfo=PATHSTR; 

    PrefixeFolder=['Hypno_' NewEXPFilename]
    HypnoDir=fullfile(paramsOldEXP.FileInfo.FilesDir,PrefixeFolder);

    mkdir(HypnoDir);

    for nHyp=1:length(NewFileInfo.HypnoFiles)
        %update the structure of the nex exp file
        NewFileInfo.HypnoFiles(nHyp).FileName=fullfile(PrefixeFolder,[NewEXPFilename '_' paramsOldEXP.FileInfo.HypnoFiles(nHyp).FileName(end-length(formatIn)-1:end)]); 
         NewFileInfo.HypnoFiles(nHyp).Dir=paramsOldEXP.FileInfo.FilesDir;
        if exist(fullfile(paramsOldEXP.FileInfo.FilesDir,NewFileInfo.HypnoFiles(nHyp).FileName))==1
            Ans=questdlg('File already exist, would you like to replace it?');
            if strcmp(Ans,'Yes')    
                %duplicate hypnofile in a new file
                copyfile(fullfile(paramsOldEXP.FileInfo.HypnoFiles(nHyp).Dir,paramsOldEXP.FileInfo.HypnoFiles(nHyp).FileName),...
                    fullfile(paramsOldEXP.FileInfo.FilesDir,NewFileInfo.HypnoFiles(nHyp).FileName));
            end
        else

             copyfile(fullfile(paramsOldEXP.FileInfo.HypnoFiles(nHyp).Dir,paramsOldEXP.FileInfo.HypnoFiles(nHyp).FileName),...
                    fullfile(paramsOldEXP.FileInfo.FilesDir,NewFileInfo.HypnoFiles(nHyp).FileName));
        end

        %update the new struct
        if iscell(s.Animal.Hypnogram.Files.File)==1
            s.Animal.Hypnogram.Files.File{nHyp}.FileName.Text=NewFileInfo.HypnoFiles(nHyp).FileName;
        elseif isstruct(s.Animal.Hypnogram.Files.File)==1
            s.Animal.Hypnogram.Files.File(nHyp).FileName.Text=NewFileInfo.HypnoFiles(nHyp).FileName;
        else
            errordlg('Error while saving the new exp file');
            break
        end

    end
    %duplicate the matfile with the new name    
    copyfile(fullfile(paramsOldEXP.FileInfo.FilesDir,[OldExpFilenameExt(1:end-3) 'mat']),...
            fullfile(paramsOldEXP.FileInfo.FilesDir,[NewEXPFilename '.mat']));



    %save the nex EXP;
    XMLCurrFileName=[NewEXPFilename '.xml'];
    struct2xml(s,fullfile(paramsOldEXP.FileInfo.FilesDir,XMLCurrFileName));
    copyfile(fullfile(paramsOldEXP.FileInfo.FilesDir,XMLCurrFileName),fullfile(paramsOldEXP.FileInfo.FilesDir,NewEXPFilenameExp),'f');
    delete(fullfile(paramsOldEXP.FileInfo.FilesDir,XMLCurrFileName)); 
    Status=1;
catch
    Status=0;
    NewFileInfo=[]
end
