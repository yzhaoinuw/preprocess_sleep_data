function ModifyHypnoFromTimeLim(NewCode)

%NewCode is the Code to add in the hypno when oscillation appears

%pick the file that contain the oscillation Lim
try
    load TempDir
    [FileName,DefaultDir,FilterIndex]=uigetfile({'*.mat'},'Pick files','multiselect','on',DefaultDir);
catch
    [FileName,DefaultDir,FilterIndex]=uigetfile({'*.mat'},'Pick files','multiselect','on','c:\');
end
if FilterIndex==1 
    load([DefaultDir FileName],'Data');
    if exist('Data')==1
        if isfield(Data,'Event')==1
           
            Oscillations=Data;

            %load the refernce expto duplicate
            FileInfo=loadEXP;


            OldExpFilenameExt=FileInfo.Filename;
            [NewEXPFilename,d]=uiputfile(fullfile(FileInfo.FilesDir,OldExpFilenameExt));

            [PATHSTR,NewEXPFilename,EXT] = fileparts([d NewEXPFilename]);
            NewFile=FileInfo;
            NewEXPFilenameExp=[NewEXPFilename '.exp'];

            formatIn = 'yyyy-mm-dd_HH-MM-SS-FFF';


            %load old exp struct
            [ s ] = Xml2struct(fullfile(FileInfo.FilesDir,OldExpFilenameExt));

            %update the fileorigin in the exp struct
            s.Animal.FileOrigin.Text=d;
            NewFile.FileInfo=PATHSTR;

            for nHyp=1:length(NewFile.HypnoFiles)
                %update the structure of the nex exp file
                NewFile.HypnoFiles(nHyp).FileName=[NewEXPFilename '_' FileInfo.HypnoFiles(nHyp).FileName(end-length(formatIn)-1:end)]; 

                %duplicate hypnofile in a new file
                copyfile(fullfile(FileInfo.HypnoFiles(nHyp).Dir,FileInfo.HypnoFiles(nHyp).FileName),...
                    fullfile(FileInfo.HypnoFiles(nHyp).Dir,NewFile.HypnoFiles(nHyp).FileName),'f');

                %update the new struct
                s.Animal.Hypnogram.Files.File{nHyp}.FileName.Text=NewFile.HypnoFiles(nHyp).FileName;


                %load HYpno
                FidH=fopen(fullfile(FileInfo.HypnoFiles(nHyp).Dir,NewFile.HypnoFiles(nHyp).FileName),'r+');
                CurrHypno=fread(FidH,'uint16');

                AbsStartLimCurrHypno=etime(datevec(FileInfo.HypnoFiles(nHyp).TStart),datevec(FileInfo.BinFiles(1).TStart));
                EndStartLimCurrHypno=AbsStartLimCurrHypno+length(CurrHypno);
                
                %modify the hypno with the newcode at specific osscilation Lim
                AllDeb=Oscillations.Event.Deb;
                AllFin=Oscillations.Event.End;

                for nOscillation=1:length(AllDeb)
                    if AllDeb(nOscillation)>=AbsStartLimCurrHypno && AllDeb(nOscillation)<=EndStartLimCurrHypno && AllFin(nOscillation)<=EndStartLimCurrHypno %if oscillation  in the current hypno
                        Deb=AllDeb(nOscillation);
                        Fin=AllFin(nOscillation);

                        StartHypnoAnalysis=floor(Deb-AbsStartLimCurrHypno);
                        HypnoCurrPart=repmat(NewCode,1,ceil(Fin-Deb));
                        WriteBinFile(1,StartHypnoAnalysis,FidH,1,uint16(HypnoCurrPart));

                    elseif AllDeb(nOscillation)>=AbsStartLimCurrHypno && AllDeb(nOscillation)<=EndStartLimCurrHypno && AllFin(nOscillation)>EndStartLimCurrHypno
                        Deb=AllDeb(nOscillation);
                        Fin=EndStartLimCurrHypno;


                        StartHypnoAnalysis=floor(Deb-AbsStartLimCurrHypno);
                        HypnoCurrPart=repmat(NewCode,1,ceil(Fin-Deb));
                        WriteBinFile(1,StartHypnoAnalysis,FidH,1,uint16(HypnoCurrPart));

                    elseif AllDeb(nOscillation)<AbsStartLimCurrHypno && AllFin(nOscillation)>=EndStartLimCurrHypno && AllFin(nOscillation)<=EndStartLimCurrHypno
                        Deb=AbsStartLimCurrHypno;
                        Fin=AllFin(nOscillation);


                        StartHypnoAnalysis=floor(Deb-AbsStartLimCurrHypno);
                        HypnoCurrPart=repmat(NewCode,1,ceil(Fin-Deb));
                        WriteBinFile(1,StartHypnoAnalysis,FidH,1,uint16(HypnoCurrPart));


                    end

                end
                fclose(FidH);


            end


            %duplicate the matfile with the new name    
            copyfile(fullfile(FileInfo.FilesDir,[OldExpFilenameExt(1:end-3) 'mat']),...
                    fullfile(FileInfo.FilesDir,[NewEXPFilename '.mat']));

            %set the current exp file with the new name
           

            %save the nex EXP;
            XMLCurrFileName=[NewEXPFilename '.xml'];
            struct2xml(s,fullfile(FileInfo.FilesDir,XMLCurrFileName));
            copyfile(fullfile(FileInfo.FilesDir,XMLCurrFileName),fullfile(FileInfo.FilesDir,NewEXPFilenameExp),'f');
            delete(fullfile(FileInfo.FilesDir,XMLCurrFileName));

        end
    end
end

'Done'





