function Artef=ArtefactLim(T0,TotDuration,ArtefThresh,Info)
%this function extract  artefact time from an ICA analysis made on slipanalysis



AbsStart=T0*60;
AbsEnd=AbsStart+TotDuration*60;
if nargin==3
    %load expfile
    Info=loadEXP;
elseif nargin==2
     %load expfile
    Info=loadEXP;
    ArtefThresh=[]; 
end



%init artefact time
ArtefTimeDeb=[];   
ArtefTimeEnd=[];

%artefact evaluation from ICA analysis
if isempty(Info.Artef)==0
    %exctract the ICA component which cointain the artefact
    IC2Remove=eval(sprintf('%s',Info.Artef.IC2Remove));

    %test if data are too large 
    nbPt2Read=round((AbsEnd-AbsStart)*Info.Fs);
    MaxPt2read=10000000 ;
    if nbPt2Read>MaxPt2read %if there is more than MaxPt2read pt to read all data are readen with a loop

        DataArtef=[];

        for nloop=1:ceil(nbPt2Read/MaxPt2read)
            ICAdata=[];

            StartLim=AbsStart+(MaxPt2read*(nloop-1))/Info.Fs;
            EndLim=AbsStart+(MaxPt2read*nloop-1)/Info.Fs;

            %check if end limit is higher than AbsEnd
            if EndLim>AbsEnd
                EndLim=AbsEnd;
            end

            %get all the data to evaluate the ICA
             BinData=ExtractContinuousData(Info.FilesDir,Info,[],StartLim,EndLim);

            % %get the chan of interest and chan 2 correct
             ChanOfInterest=eval(Info.Artef.ChanOfInterest); 

             %compute ICA on current data
             ICAdata=Info.Artef.ICAMat(:,:,1)*BinData(ChanOfInterest,:);

             DataArtef=[DataArtef ICAdata(find(ChanOfInterest==IC2Remove),:)];


             clear BinData ICAdata ;

        end


    else %if there is less than MaxPt2read pt to read all data without any loop

        %get all the data to evaluate the ICA
         BinData=ExtractContinuousData(Info.FilesDir,Info,[],AbsStart,AbsEnd);

        % %get the chan of interest and chan 2 correct
         ChanOfInterest=eval(Info.Artef.ChanOfInterest); 

         %compute ICA on current data
         ICAdata=Info.Artef.ICAMat(:,:,1)*BinData(ChanOfInterest,:);
         DataArtef=ICAdata(find(ChanOfInterest==IC2Remove),:);
         clear BinData ICAdata;
    end



%              %remove mean
    DataArtef=abs(DataArtef-nanmean(DataArtef));


    %detect limit time 

    %extraction du max par fenetre glissante de 10s

    WinArtef=round(10*Info.Fs);
    x = vanherk(DataArtef,WinArtef,'max','same');
    
    %plot artefact detection
    hArtef=figure;
    LinePlotReducer(@plot,DataArtef);hold on;LinePlotReducer(@plot,x,'r');
    
    
    if isempty(ArtefThresh)==1
        ArtefThreshStr=inputdlg('ArtefThresh=');
        ArtefThresh=str2double(ArtefThreshStr{1});
    end

    idxArtef=abs(x)>ArtefThresh;

    
    LinePlotReducer(@plot,[1 length(x)],[ArtefThresh ArtefThresh],'k');
    grid on;
    
    set(hArtef,'units','normalized');
    set(hArtef,'units','normalized','position',[0 0 1 1]);
    
    
    

    %close(gcf);
    Fulltime=linspace(AbsStart,AbsEnd,length(idxArtef)+1);

    %extract artefact start and end times
    ArtefTimeDeb=Fulltime(diff([0 idxArtef 0])==1);
    ArtefTimeEnd=Fulltime(diff([0 idxArtef 0])==-1);

    Artef.Thresh=ArtefThresh;

    %save artefact
    Artef.ArtefTimeDeb=ArtefTimeDeb;
    Artef.ArtefTimeEnd=ArtefTimeEnd;
    
    ArtefSetting=Info.Artef;
    
    %save mat file
    MatFileName=sprintf('%s%s ICA N°%d Thresh %d %06.1f min-%06.1f min.mat',Info.FilesDir,Info.ExpFileName(1:end-4),IC2Remove,ArtefThresh,T0,T0+TotDuration);
    save(MatFileName,'Artef','ArtefSetting');
    
    %save figure
    
    ArtefFigFileName=sprintf('%s%s ICA N°%d Thresh %d %06.1f min-%06.1f min.fig',Info.FilesDir,Info.ExpFileName(1:end-4),IC2Remove,ArtefThresh,T0,T0+TotDuration);
    %saveas(hArtef,ArtefFigFileName,'fig');
    
    hgsave(hArtef, ArtefFigFileName, '-v7.3')
    
    
    
    ArtefTifFileName=sprintf('%s%s ICA N°%d Thresh %d %06.1f min-%06.1f min.tif',Info.FilesDir,Info.ExpFileName(1:end-4),IC2Remove,ArtefThresh,T0,T0+TotDuration);
    
   
    set(hArtef,'PaperPositionMode','Auto');
   
    print(hArtef,ArtefTifFileName,'-dtiff','-r100');
    
    %close(hArtef);
    


    
else

    'No Artefact Info in the Mat file!'
    Artef=[];


end 

