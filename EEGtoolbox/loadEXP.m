function Info=loadEXP(varargin) %% load exp file

%First argument should be the full filename and the second 'Yes' or 'No'
%depending if the AcquisitionRangeMax (in the exp file) is half the real data acquisition range

%varargin{1} is the full filename
%varagin{2} is 'Yes' or 'No' answer to the question about the acquisition range max
  
    if isempty(varargin) || (isempty(varargin)==0 && isempty(varargin{1}))

        try
            load TempDir
            [FileName,DefaultDir,FilterIndex]=uigetfile({'*.exp'},'Pick files','multiselect','off',DefaultDir);
        catch
            [FileName,DefaultDir,FilterIndex]=uigetfile({'*.exp'},'Pick files','multiselect','off','c:\');
        end

    elseif isempty(varargin{1})==0
        filename=varargin{1};
        [DefaultDir,FileName,ext]=fileparts(filename);
        DefaultDir=[DefaultDir '\'];
        FileName=[FileName ext];
        FilterIndex=1;

    end
    
    if strcmp(DefaultDir(end),'\')==1
        DefaultDir=DefaultDir(1:end-1);
    end


    if nargin<2 || (nargin==2 && isempty(varargin{2}))
        MaxRangeType=questdlg('The AcquisitionRangeMax is the half the whole acquisition range?')  ;  
    else
        MaxRangeType=varargin{2};
    end
     if FilterIndex==1
 
        try

            save TempDir DefaultDir;
            %save the parameters into Info

            Info.FilesDir=DefaultDir;
            Info.ExpFileName=FileName;
            

            
            %create variable with all filename
            [FileDir,FileWithoutExt,ext] = fileparts(FileName);
            
            allFileName = fullfile(DefaultDir,[FileWithoutExt ext]);
    

            %Read all exp file
            [ s ] = Xml2struct(allFileName);
            Channels=s.Animal.Acquisition.Channels.Channel;
            ChLabel={};
            Gain=[];
            Offset=[];
            Range=[];
            
            NbActivChan=0;

            for nchan=1:length(Channels)
                try
                    CurrChS=Channels{nchan};
                catch
                    CurrChS=Channels(nchan);
                end
                
                if strcmp(CurrChS.Enable.Text,'true')==1 %Update if channel is enable
                    NbActivChan=NbActivChan+1;
                    ChLabel=[ChLabel CurrChS.Name.Text];
                    Offset=[Offset str2double(CurrChS.Offset.Text)];
                    Gain=[Gain str2double(CurrChS.Gain.Text)];
                    if strcmp(MaxRangeType,'Yes')==1
                         %the acquisitionrange max is half the whole range
                        Range=[Range str2double(CurrChS.AcquisitionRangeMax.Text).*2];
                    else
                        Range=[Range str2double(CurrChS.AcquisitionRangeMax.Text)];
                    end
                end
%                 Range=[Range str2double(CurrChS.AcquisitionRangeMax.Text)];
            end


            % Paramètres d'acquisition
            Info.NbRecChan=NbActivChan;
            Info.Name=s.Animal.Name.Text;
           % Info.NbRecChan = str2double(s.Animal.Acquisition.NbChan.Text);
            Info.ChLabel = ChLabel;
            Info.Fs = str2double(s.Animal.Acquisition.SamplingRate.Text);
            Info.Gain = Gain;
            Info.Range = Range;
            Info.Offset =Offset; % le point 0 du 16 bit
            try
                Info.Precision = s.Animal.Acquisition.AcquisitionType.Text;
                if strcmpi(Info.Precision,'uint16')==1
                   Info.Precision= 'uint16';
                elseif strcmpi(Info.Precision,'int16')==1
                    Info.Precision= 'int16';
                end
            catch
                Info.Precision = 'uint16';
            end

            Info.Filename = [FileWithoutExt ext];

            
            %init Artef and State
            Info.State=[];
            Info.Artef=[];
            
            %try to load the state from the exp file
            try
                
                ExpStateList=s.Animal.Hypnogram.Statuses.Status;
                for nstate=1:length(ExpStateList)
                    Info.State(nstate).Label=ExpStateList{nstate}.Label.Text;
                    Info.State(nstate).Code=str2double(ExpStateList{nstate}.Level.Text);
                    if length(ExpStateList{nstate}.Key.Text)==1 && str2double(ExpStateList{nstate}.Key.Text)<=9
                        Info.State(nstate).Key=['numpad' ExpStateList{nstate}.Key.Text];
                    else
                        Info.State(nstate).Key=ExpStateList{nstate}.Key.Text;
                    end
                    Info.State(nstate).Color=Dec2ARGB(str2double(ExpStateList{nstate}.Color.Text))./255;
  
                end
          
            end
            
            
            %try to importe the states value and color from the mat file
            %mat file has the priority to the exp file
            try 
                AnimalMatFileName=[allFileName(1:end-3) 'mat'];
                load(AnimalMatFileName);
                
                try
                    Info.State=PStates;
                catch
                    Info.State=Animal.PStates;
                end
                 
                try %try to load ICA detection
                    Info.Artef=[];
                    Info.Artef=PArtef;
 
                catch
                    Info.Artef=PStates;
                    
                end

            end
           
                

            %bin files
            BinFiles=s.Animal.Acquisition.Files.File;

            if length(BinFiles)==1 %if only one bin file there is no binfile idx is not a cell
                BinFiles={BinFiles}; 
            end
            
            

            
            for nbin=1:length(BinFiles)
                %get only the filename without the parent folder if exist
                %in the exp file
                
                

                
                formatOut = 'yyyy-mm-ddTHH:MM:SS.FFF';
                try
                    Info.BinFiles(nbin).TStart=datenum(BinFiles{nbin}.TStart.Text,formatOut);%Start in relative time
                catch
                    Info.BinFiles(nbin).TStart=datenum([BinFiles{nbin}.TStart.Text '.000'],formatOut);%Start in relative time
                end
                
                %create sub dir filed that could containt the subfolder if bin are stored in a subfolder
                [~,filename,ext]=fileparts(BinFiles{nbin}.FileName.Text);
                    Info.BinFiles(nbin).FileName=[filename ext];
                
                if exist(fullfile(DefaultDir,BinFiles{nbin}.FileName.Text))==2;
                    [Info.BinFiles(nbin).Dir,~]=fileparts(fullfile(DefaultDir,BinFiles{nbin}.FileName.Text));
                    
                else
                    
                    
                    try
                        Out=recursdir(DefaultDir,Info.BinFiles(nbin).FileName);
                        Info.BinFiles(nbin).Dir=Out.Dir{1};                
                    catch
                        Info.BinFiles(nbin).Dir='';
                    end
                end
                
                try
                    Info.BinFiles(nbin).Duration=GetBinDuration(fullfile(Info.BinFiles(nbin).Dir,Info.BinFiles(nbin).FileName),Info.NbRecChan)/Info.Fs;
                catch
                    'Empty bin file or GetBinDuration.m is missing'
                    Info.BinFiles(nbin).Duration=0;
                end
            end
            
            
%%%%%%%%%%%%%%%%%%%%%
%Video and actimetry files
            
%             
%              if isfield(s.Animal,'Videos')==1 && isfield(s.Animal.Videos,'Video')==1
%                 Videos=s.Animal.Videos.Video;
%                 NbVideos=length(Videos);
%                 if NbVideos==1
%                     Videos={Videos};
%                 end
% 
%                 for Nv=1:NbVideos %for each videos
% %                     VideosFiles=[];
%                     try   
% 
% %                         Info.VideosFiles(Nv).FrameRate=str2double(Videos{Nv}.SamplingRate.Text);
% % 
% % 
%                         VideosFiles=Videos{Nv}.Files.File;
%                         if length(VideosFiles)==1
%                             VideosFiles={VideosFiles};
%                         end
% 
%                         try
%                             if isfield(Videos{Nv},'Actimetry') && isfield(Videos{Nv}.Actimetry,'Files')==1 && isfield(Videos{Nv}.Actimetry.Files,'File')==1
%                                 ActiFiles=Videos{Nv}.Actimetry.Files.File;
%                                 
%                        
% 
% 
%                                  if length(ActiFiles)==1
%                                     ActiFiles={ActiFiles};
%                                  end
%                                  
%                                 %remove empty act file                                          
%                                 Idok=false(length(ActiFiles),1);
%                                 for nActi=1:length(ActiFiles);
%                                     ActifileName=[DefaultDir Videos{Nv}.Actimetry.Files.File{nActi}.FileName.Text];
%                                     fid=fopen(ActifileName);
%                                     NbSample=length(fread(fid,[2,inf],'int32'));
%                                     fclose(fid);
%                                     
%                                     if NbSample>0
%                                         Idok(nActi)=true;
%                                 
%                                     end
%                                 end
%                                 ActiFiles=ActiFiles(Idok);
%                                  
% 
%                                  Info.ActiFiles(Nv).SamplingRate=str2double(Videos{Nv}.Actimetry.SamplingRate.Text);
%                                  Info.ActiFiles(Nv).Name=Videos{Nv}.Actimetry.Name.Text;
%                                  
%                                  
%                                   %check if ActiFiles is empty
%                                 for nacti=1:length(ActiFiles)
%                                     try
%                                         NbSample=GetBinDuration(fullfile(DefaultDir,ActiFiles{nacti}.FileName.Text),2);
%                                         if NbSample==0
%                                             ActiFiles{nacti}=[];
% 
%                                         end
%                                     catch
%                                          ActiFiles{nacti}=[];
%                                     end
% 
% 
%                                 end
% 
% 
%                             else
%                                 ActiFiles=[];
%                             end
%                         catch ME1
% 
%                             ME1.identifier
%                             ME1.message
%                             for i=1:length(ME1.stack)
%                                 ME1.stack(i)
%                             end
%                         end
% 
% 
% 
%                          for nvideo=1:length(ActiFiles)
% 
%                             %try to open the disp actimetry file
%                             try  
%                                 if isempty(ActiFiles{nvideo})==0
%                                      Info.ActiFiles(Nv).Files(nvideo).FileName=ActiFiles{nvideo}.FileName.Text;
%                                      formatOut = 'yyyy-mm-ddTHH:MM:SS.FFF';
%                                     try
%                                         Info.ActiFiles(Nv).Files(nvideo).TStart=datenum(ActiFiles{nvideo}.TStart.Text,formatOut);%Start in relative time
%                                     catch
%                                         Info.ActiFiles(Nv).Files(nvideo).TStart=datenum([ActiFiles{nvideo}.TStart.Text '.000'],formatOut);%Start in relative time
%                                     end
%                                     Info.ActiFiles(Nv).Files(nvideo).Duration=str2double(ActiFiles{nvideo}.Duration.Text)/1000; %duration in s
%                                 end
%                             end
% 
%                          end
% 
% 
% 
%                          %reorder the file if not in alphabetical order               
%                         % [ListName,idx]=sort({Info.ActiFiles(Nv).Files(:).FileName});
%                          ActiName={Info.ActiFiles(Nv).Files(:).FileName};
%                          ActiStart=[Info.ActiFiles(Nv).Files(:).TStart];
%                          ActiDuration=[Info.ActiFiles(Nv).Files(:).Duration];
%                         [~,idx]=sort(ActiStart);
% 
%                          for nidx=1:length(idx)
% 
%                              Info.ActiFiles(Nv).Files(nidx).FileName=ActiName{idx(nidx)};
%                              Info.ActiFiles(Nv).Files(nidx).TStart=ActiStart(idx(nidx));
%                              Info.ActiFiles(Nv).Files(nidx).Duration=ActiDuration(idx(nidx));
% 
%                          end
% 
%                          
%                       catch ME1
% 
%                       ME1.identifier
%                       ME1.message
%                       for i=1:length(ME1.stack)
%                           ME1.stack(i)
%                       end   
% 
% 
%                     end
%                 
%                         
%                 end     
%                         
%                         
%                         
%                         
%                 %zoom(handles.DataActimetry,'on');
%             else
%                 % Info.VideosFiles=[]; 
%                  Info.ActiFiles=[];
% 
%             end
% 
% 
%             
     


%read video
            try
                Info.VideosFiles=[];
                Info.ActiFiles=[];
                if isfield(s.Animal,'Videos')==1 && isfield(s.Animal.Videos,'Video')==1
                    Videos=s.Animal.Videos.Video;
                    NbVideos=length(Videos);
                    if NbVideos==1
                        Videos={Videos};
                    end
                    
                    NbActimetry=1;
                    
                    for Nv=1:NbVideos %for each videos

                        try

                                Info.VideosFiles(Nv).FrameRate=str2double(Videos{Nv}.SamplingRate.Text);
                                Info.VideosFiles(Nv).Name=Videos{Nv}.Name.Text;

                                VideosFiles=Videos{Nv}.Files.File;
                                if length(VideosFiles)==1
                                    VideosFiles={VideosFiles};
                                end
                                

%                                 %extraction des info pour chaque fichier d'actimetrie de la video courante
                                try
                                    if isfield(Videos{Nv},'Actimetry') && isfield(Videos{Nv}.Actimetry,'Files')==1 && isfield(Videos{Nv}.Actimetry.Files,'File')==1
                                        ActiFiles=Videos{Nv}.Actimetry.Files.File;

                                         if length(ActiFiles)==1
                                            ActiFiles={ActiFiles};
                                         end



                                      %extract actifile info                                         
                                        %Idok=false(length(ActiFiles),1);
                                        NbActiFile(Nv)=0;
                                        for nActi=1:length(ActiFiles);

                                            %create sub dir filed that could containt the subfolder if bin are stored in a subfolder
                                            [~,filename,ext]=fileparts( ActiFiles{nActi}.FileName.Text);
                                             filename=[filename ext];
                                            
                                            if exist(fullfile(DefaultDir,ActiFiles{nActi}.FileName.Text))==2;
                                                [ActiFilesDir,~]=fileparts(fullfile(DefaultDir,ActiFiles{nActi}.FileName.Text));
                                       
                                            else

                                                
                                                try
                                                   
                                                    Out=recursdir(DefaultDir,filename);
                                                    ActiFilesDir=Out.Dir{1};                
                                                catch
                                                    ActiFilesDir='';
                                                end
                                            end
%                                             
%                                             
%                                             
%                                             [~,filename,ext]=fileparts( ActiFiles{nActi}.FileName.Text);
%                                             ActiFiles{nActi}.FileName.Text=[filename ext];
%                                           
%                                             
%                                             Out(nActi)=recursdir(DefaultDir,ActiFiles{nActi}.FileName.Text);
                                            ActifileName=fullfile(ActiFilesDir,filename);
                                            fid=fopen(ActifileName);
                                            
                                            if fid>0
                                                NbSample=length(fread(fid,[2,inf],'int32'));
                                                fclose(fid);
                                            else
                                                NbSample=0;
                                            end 

                                            if NbSample>0
                                                
                                                NbActiFile(Nv)=NbActiFile(Nv)+1;
                                                %Idok(nActi)=true;
                                                 Info.ActiFiles(Nv).Files(NbActiFile(Nv)).FileName=filename;
                                                 %Out=recursdir(DefaultDir,ActiFiles{nActi}.FileName.Text);
                                                 Info.ActiFiles(Nv).Files(NbActiFile(Nv)).Dir=ActiFilesDir;
                                                 formatOut = 'yyyy-mm-ddTHH:MM:SS.FFF';
                                                try
                                                    Info.ActiFiles(Nv).Files(NbActiFile(Nv)).TStart=datenum(ActiFiles{nActi}.TStart.Text,formatOut);%Start in relative time
                                                catch
                                                    Info.ActiFiles(Nv).Files(NbActiFile(Nv)).TStart=datenum([ActiFiles{nActi}.TStart.Text '.000'],formatOut);%Start in relative time
                                                end
                                                Info.ActiFiles(Nv).Files(NbActiFile(Nv)).Duration=str2double(ActiFiles{nActi}.Duration.Text)/1000; %duration in s
%                                             else
%                                                 Info.ActiFiles(Nv).Files(nActi).Name='';
%                                                 Info.ActiFiles(Nv).Files(nActi).Duration='';
%                                                 Info.ActiFiles(Nv).Files(nActi).TStart=[];
                                            end
                                        end
                                        
                                        
                                        if NbActiFile(Nv)>=0 %if at least 1 file exst
                                           
                                            
                                           % ActiFiles=ActiFiles(Idok);
                                            %Out=Out(Idok);
                                            Info.ActiFiles(NbActimetry).SamplingRate=str2double(Videos{Nv}.Actimetry.SamplingRate.Text);
                                            Info.ActiFiles(NbActimetry).Name=Videos{Nv}.Name.Text;
                                            
                                             %reorder the file if not in alphabetical order               
                                            % [~,idx]=sort({Info.VideosFiles(Nv).Files(:).FileName});
                                             ActName={Info.ActiFiles(NbActimetry).Files(:).FileName};
                                             ActDir={Info.ActiFiles(NbActimetry).Files(:).Dir};
                                             ActStart=[Info.ActiFiles(NbActimetry).Files(:).TStart];
                                             ActDuration=[Info.ActiFiles(NbActimetry).Files(:).Duration];

                                             [~,idx]=sort(ActStart);


                                             for nidx=1:length(idx)

                                                 Info.ActiFiles(NbActimetry).Files(nidx).FileName=ActName{idx(nidx)};
                                                 Info.ActiFiles(NbActimetry).Files(nidx).Dir=ActDir{idx(nidx)};
                                                 Info.ActiFiles(NbActimetry).Files(nidx).TStart=ActStart(idx(nidx));
                                                 Info.ActiFiles(NbActimetry).Files(nidx).Duration=ActDuration(idx(nidx));

                                             end


                                        end
                                        

                                        
                                        
                                        
                                        
                                        
%                                         
%                                         %  check if ActiFiles is empty
%                                         for nacti=1:length(ActiFiles)
% 
%                                             try
%                                                 NbSample=GetBinDuration(fullfile(Out(nacti).Dir{1},ActiFiles{nacti}.FileName.Text),2);
%                                                 if NbSample==0
%                                                     ActiFiles{nacti}=[];
% 
%                                                 end
%                          
%                                             catch
%                                                  ActiFiles{nacti}=[];
%                                             end
% 
% 
%                                         end
                                        NbActimetry=NbActimetry+1;
                                    end
                                catch ME1

                                    ME1.identifier
                                    ME1.message
                                    for i=1:length(ME1.stack)
                                        ME1.stack(i)
                                    end
                                end


                                 nvidok=0;
                                 for nvideo=1:length(VideosFiles)

%                                      [~,filename,ext]=fileparts(VideosFiles{nvideo}.FileName.Text);
%                                      Info.VideosFiles(Nv).Files(nvideo).FileName=[filename ext];
%                                    
%                                     
%                                     Out=recursdir(DefaultDir,Info.VideosFiles(Nv).Files(nvideo).FileName);         
%                                     Info.VideosFiles(Nv).Files(nvideo).Dir=Out.Dir{1};

                                    %create sub dir filed that could containt the subfolder if bin are stored in a subfolder
                                    [~,filename,ext]=fileparts(VideosFiles{nvideo}.FileName.Text);
                                    if nvideo==1 || (nvideo>1 && sum(strcmp({Info.VideosFiles(Nv).Files(:).FileName},[filename ext]))==0)
                                            nvidok=nvidok+1;
                                            Info.VideosFiles(Nv).Files(nvideo).FileName=[filename ext];
                                       
                                   
                                    
                                        if exist(fullfile(DefaultDir,VideosFiles{nvideo}.FileName.Text))==2;
                                            [Info.VideosFiles(Nv).Files(nvideo).Dir,~]=fileparts(fullfile(DefaultDir,VideosFiles{nvideo}.FileName.Text));
                                        else


                                            try
                                                Out=recursdir(DefaultDir,Info.VideosFiles(Nv).Files(nvideo).FileName);
                                                Info.VideosFiles(Nv).Files(nvideo).Dir=Out.Dir{1};                
                                            catch
                                                Info.VideosFiles(Nv).Files(nvideo).Dir='';
                                            end
                                        end






                                        formatOut = 'yyyy-mm-ddTHH:MM:SS.FFF';
                                        try
                                            Info.VideosFiles(Nv).Files(nvideo).TStart=datenum(VideosFiles{nvideo}.TStart.Text,formatOut);%Start in relative time
                                        catch
                                            Info.VideosFiles(Nv).Files(nvideo).TStart=datenum([VideosFiles{nvideo}.TStart.Text '.000'],formatOut);%Start in relative time
                                        end
                                        Info.VideosFiles(Nv).Files(nvideo).Duration=str2double(VideosFiles{nvideo}.Duration.Text)/1000; %duration in s


                                        %get the frame rate from file header

                                        try
                                            VidInfo=extractAviInfo(fullfile(Info.VideosFiles(Nv).Files(nvideo).Dir,Info.VideosFiles(Nv).Files(nvideo).FileName));
                                            Info.VideosFiles(Nv).Files(nvideo).FrameRate=VidInfo.Fps;
                                            Info.VideosFiles(Nv).Files(nvideo).TimeStamp=[];
                                            %check if there is a file containing th etimestamp
                                            TimeStampFileName=fullfile(Info.VideosFiles(Nv).Files(nvideo).Dir,[Info.VideosFiles(Nv).Files(nvideo).FileName(1:end-4) '_Timers.txt']);
                                            matTimeStampFileName=fullfile(Info.VideosFiles(Nv).Files(nvideo).Dir,[Info.VideosFiles(Nv).Files(nvideo).FileName(1:end-4) '_Timers.mat']);

                                            if exist(TimeStampFileName) && exist(matTimeStampFileName)==0


                                              [ImTimeStamp,ImNum,ImErr,meanfps,totalfps]=extractTimstampFromTxtFile(TimeStampFileName);                                          Info.VideosFiles(Nv).Files(nvideo).FrameRate=totalfps; 



                                              Info.VideosFiles(Nv).Files(nvideo).TimeStamp=ImTimeStamp;
                                              Info.VideosFiles(Nv).Files(nvideo).FrameRate=totalfps;
                                              save(matTimeStampFileName,'ImTimeStamp','ImNum','ImErr','meanfps','totalfps');
                                              sprintf('%s, loaded',[Info.VideosFiles(Nv).Files(nvideo).FileName(1:end-4) '_Timers.txt'])

                                            elseif exist(TimeStampFileName) && exist(matTimeStampFileName)

                                              load(matTimeStampFileName,'ImTimeStamp','totalfps'); ;    
                                              Info.VideosFiles(Nv).Files(nvideo).TimeStamp=ImTimeStamp;
                                              Info.VideosFiles(Nv).Files(nvideo).FrameRate=totalfps;
                                              sprintf('%s, loaded',[Info.VideosFiles(Nv).Files(nvideo).FileName(1:end-4) '_Timers.txt'])
                                            end

                                        catch
                                            'impossible to read the frame rate into the header of the avi file!'
                                            if isfield(VideosFiles{nvideo},'SamplingRate') %if a specific sampling rate from the video exist
                                                Info.VideosFiles(Nv).Files(nvideo).FrameRate=str2double(VideosFiles{nvideo}.SamplingRate.Text);
                                                Info.VideosFiles(Nv).Files(nvideo).TimeStamp=[];
                                            else
                                                Info.VideosFiles(Nv).Files(nvideo).FrameRate=Info.VideosFiles(Nv).FrameRate;
                                                Info.VideosFiles(Nv).Files(nvideo).TimeStamp=[];
                                            end
                                        end

                                    end

                                 end



                                 %reorder the file if not in alphabetical order               
                                % [~,idx]=sort({Info.VideosFiles(Nv).Files(:).FileName});
                                 VidName={Info.VideosFiles(Nv).Files(:).FileName};
                                 VidDir={Info.VideosFiles(Nv).Files(:).Dir};
                                 VidStart=[Info.VideosFiles(Nv).Files(:).TStart];
                                 VidDuration=[Info.VideosFiles(Nv).Files(:).Duration];

                                 [~,idx]=sort(VidStart);


                                 for nidx=1:length(idx)

                                     Info.VideosFiles(Nv).Files(nidx).FileName=VidName{idx(nidx)};
                                     Info.VideosFiles(Nv).Files(nidx).Dir=VidDir{idx(nidx)};
                                     Info.VideosFiles(Nv).Files(nidx).TStart=VidStart(idx(nidx));
                                     Info.VideosFiles(Nv).Files(nidx).Duration=VidDuration(idx(nidx));

                                 end

                          catch ME1

                        ME1.identifier
                        ME1.message
                        for i=1:length(ME1.stack)
                            ME1.stack(i)
                        end   


                        end
                    end



                else
                     Info.VideosFiles=[]; 
                     Info.ActiFiles=[];

                end

            catch ME1

                ME1.identifier
                ME1.message
                for i=1:length(ME1.stack)
                    ME1.stack(i)
                end
                  Info.VideosFiles=[];   
                  Info.ActiFiles=[];

             end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %Hypno files
            try
                HypnoFiles=s.Animal.Hypnogram.Files.File;
                if length(HypnoFiles)==1 %if only one bin file there is no binfile idx is not a cell
                    HypnoFiles={HypnoFiles};
                end


                for nhyp=1:length(HypnoFiles)
%                     [~,filename,ext]=fileparts(HypnoFiles{nhyp}.FileName.Text);
%                     Info.HypnoFiles(nhyp).FileName=[filename ext];
                    try
                        formatOut = 'yyyy-mm-ddTHH:MM:SS';
                        Info.HypnoFiles(nhyp).TStart=datenum(HypnoFiles{nhyp}.TStart.Text,formatOut);%Start in relative time
                    catch  %import data in num if not converted with formatOut
                        Info.HypnoFiles(nhyp).TStart=str2double(HypnoFiles{nhyp}.TStart.Text);%Start in relative time
                    end
                    
                    
                    %create sub dir filed that could containt the subfolder if bin are stored in a subfolder
                    [~,filename,ext]=fileparts(HypnoFiles{nhyp}.FileName.Text);
                        Info.HypnoFiles(nhyp).FileName=[filename ext];
                    
                    if exist(fullfile(DefaultDir,HypnoFiles{nhyp}.FileName.Text))==2;
                        [Info.HypnoFiles(nhyp).Dir,~]=fileparts(fullfile(DefaultDir,HypnoFiles{nhyp}.FileName.Text));
                    else

                        
                        try
                            Out=recursdir(DefaultDir,Info.HypnoFiles(nhyp).FileName);
                            Info.HypnoFiles(nhyp).Dir=Out.Dir{1};                
                        catch
                            Info.HypnoFiles(nhyp).Dir='';
                        end
                    end
                    
                    
                    
                    
                    try



                        %Out=recursdir(Info.FilesDir,Info.HypnoFiles(nhyp).FileName);

                        NumberOfSamples=GetBinDuration(fullfile(Info.HypnoFiles(nhyp).Dir,Info.HypnoFiles(nhyp).FileName),1);
                        Info.HypnoFiles(nhyp).Duration=NumberOfSamples;
                       % Info.HypnoFiles(nhyp).Dir=Out.Dir{1};
                    catch ME1
                        
                          Info=[];
                          ME1.identifier
                          ME1.message
                           for i=1:length(ME1.stack)
                               ME1.stack(i)
                           end
                        
                        'Hypno file or GetBinDuration.m is missing'
                        Info.HypnoFiles(nhyp).Duration=0;
                        Info.HypnoFiles(nhyp).Dir='';
                    end

                end
                
                
                 %reorder the hypnofile if not in alphabetical order 
                try
                    HypName={Info.HypnoFiles(:).FileName};
                    HypStart=[Info.HypnoFiles(:).TStart];
                    HypDuration=[Info.HypnoFiles(:).Duration];



                    [~,idx]=sort(HypStart);

                    for nidx=1:length(idx)

                        Info.HypnoFiles(nidx).FileName=HypName{idx(nidx)};
                        Info.HypnoFiles(nidx).TStart=HypStart(idx(nidx));
                        Info.HypnoFiles(nidx).Duration=HypDuration(idx(nidx));

                    end
                    
                    %get bouts duration information
                    Info.BoutsDuration=1/str2double(s.Animal.Hypnogram.SamplingRate.Text);
           

                catch
                    'error when reordering hypno files'
                end
                
                
                
            catch
                'error when importing hypno'
            end


            %reorder the file if not in alphabetical order 

           
            try
                BinName={Info.BinFiles(:).FileName};
                BinStart=[Info.BinFiles(:).TStart];
                BinDuration=[Info.BinFiles(:).Duration];
               

                [~,idx]=sort(BinStart);

                for nidx=1:length(idx)

                    Info.BinFiles(nidx).FileName=BinName{idx(nidx)};
                    Info.BinFiles(nidx).TStart=BinStart(idx(nidx));
                    Info.BinFiles(nidx).Duration=BinDuration(idx(nidx));


                end
                
                
                
            catch
                'error when reordering bin files'
            end

            'Exp file loaded'
        catch ME1
            Info=[];
            ME1.identifier
            ME1.message
            for i=1:length(ME1.stack)
                ME1.stack(i)
            end
            warndlg('exp file not loaded');
            return;
        end
     end
end

function Out=recursdir(baseDir,filename)
% OUT = RECURSDIR(BASEDIRECTORY,SEARCHEXPRESSION)
% A recursive search to find files that match the search expression
%
%  filename
%  tic
% [status,list]=system(['dir /S ' sprintf('"%s\\*%s"',baseDir,filename)]);
% 
% %convert in line
% endofline=strfind(list, sprintf('\n'));
% s=1;
% Out={};
% nbline=0;
% for n=1:length(endofline)
%    Currl=list(s:endofline(n)-1);
%    if isempty(Currl)==0
%        
%        
%        %search for the directory
%        pos=strfind(Currl,baseDir);
%        if isempty(pos)==0 %if the directory exist (mean that the file has been found)
%            nbline=nbline+1;
%           %lengthDir=length(baseDir);
%             Out.Dir{nbline}=Currl(pos:end);
%             Out.files{nbline}=filename;
%        end
%    end
%     s=endofline(n)+1;
% end



Out.files = {};
Out.Dir= {};
% DDir=dir(baseDir);
% DDir=DDir([DDir(:).isdir]==1);
%find directory
%filter the extention file

% if nargin==3
%     dstr = dir(fullfile(baseDir,['*.' ext]));
%     dstr=cat(1,DDir,dstr);
% 
% else
    dstr = dir(baseDir);%search current directory and put results in structure
% end

%sort the dir at the end
IdDir=[dstr(:).isdir];
[~,idx]=sort(IdDir);


    
    
dstr=dstr(idx);
for II = 1:length(dstr)
    if ~dstr(II).isdir && strcmp(dstr(II).name,filename)
    %look for a match that isn't a directory 
        Out.files{length(Out.files)+1} = dstr(II).name;
        Out.Dir{length(Out.Dir)+1}=baseDir;
        break;
    elseif dstr(II).isdir && ~strcmp(dstr(II).name,'.') && ~strcmp(dstr(II).name,'..') 
    %if it is a directory(and not current or up a level), search in that
        pname = fullfile(baseDir,dstr(II).name);
        OutfilesTemp=recursdir(pname,filename);
        if ~isempty(OutfilesTemp.files)
        %if recursive search is fruitful, add it to the current list
            Out.files((length(Out.files)+1):(length(Out.files)+length(OutfilesTemp.files))) = OutfilesTemp.files;
            Out.Dir((length(Out.Dir)+1):(length(Out.Dir)+length(OutfilesTemp.Dir))) = OutfilesTemp.Dir;
            break;
        end
    end
end

% 
if isempty(Out.Dir)
    Out.Dir={''};
end
% 
if length(Out.files)>1
    warndlg('Multiple files found');
end

end