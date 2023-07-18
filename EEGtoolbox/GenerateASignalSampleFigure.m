function GenerateASignalSampleFigure(Pref,TimeStarts, Duration, Channels,XScale,YScale,HP,LP,PixH, Filename,YZero)
%Pref prefix
%TimeStarts is in sec the start of each signal, TimeStarts contain as many
%value as the number of signal to display
%Duration in second, the duration of the signal to display
%Channels channels number to display could be a differential cahnnels is a
%[2*nchan] matrix

%XScale and YScale are the length of the reference line in s for X and V or g or °C for Y
%the Yline will be 5 time less than the full axis amplitude
%HP and LP could eb enpty but this input containt the cutoff frequency
%respectively for each chann, use Nan if no filter
%PixH by default is 640 height of the output figue in pix
%Filename is the full director of the expfile to use could be empty 
%YZero the middle value on the y axe if Nan mean use a DCremove

 if nargin<=8
     PixH=640;
     Info=loadEXP([],'No'); %load the exp file to process
     YZero=NaN(1,length(Channels));
 elseif nargin<=9
     Info=loadEXP([],'No');
     YZero=NaN(1,length(Channels));
 elseif nargin<=10
     YZero=NaN(1,length(Channels));
 else
     Info=loadEXP(Filename,'No');
 end
     
CurrDir=Info.FilesDir;

nHaxe=size(Channels,2);
if size(Channels,1)==1
    Channels=[Channels;Channels];
end
nWaxe=length(TimeStarts);

Fig=figure('name','Raw Signal Sample');

PageRatio=21/29.7;
HighPixel=PixH;
set(Fig,'units','pixel','position',[20 20 HighPixel/PageRatio HighPixel ],'Color',[1 1 1]);
set(Fig,'units','normalized');

OffsetL=0.05;
OffsetR=0.01;
WSpace=0.01;
OffsetU=0.01;
OffsetD=0.05;
HSpace=0.01;

if length(YScale)==1
    YScale=repmat(YScale,1,length(Channels));
end


Width=(1-OffsetL-OffsetR+WSpace)/nWaxe-WSpace;
High=(1-OffsetD-OffsetU+HSpace)/nHaxe-HSpace;

Channels=fliplr(Channels);
YScale=fliplr(YScale);
HP=fliplr(HP);
LP=fliplr(LP);
YZero=fliplr(YZero);


for nChan=1:nHaxe
    if isempty(HP)==0 && isnan(HP(nChan))==0
        dHp = fdesign.highpass('n,fc',10,HP(nChan),Info.Fs);
        DHp = design(dHp,'butter');
        
    else
        DHp=[];
    end
    
    if isempty(LP)==0 && isnan(LP(nChan))==0
        dLp = fdesign.lowpass('n,fc',10,LP(nChan),Info.Fs);
        DLp = design(dLp,'butter');
    else
       DLp=[];
    end
    
    for nstart=1:nWaxe
        
        RelStart=TimeStarts(nstart);
        RelEnd=RelStart+Duration;

        
        %Get The mode of the state
        [ModeState,~]=GetStateBetweenLim(Info,RelStart,Duration);
        CurrColIdx=find([Info.State.Code]==ModeState);
        CurrCol=Info.State(CurrColIdx).Color;
        
    
        X=OffsetL+(nstart-1)*(Width+WSpace);
        Y=OffsetU+(nChan-1)*(High+HSpace);
        H=High;
        W=Width;
        h(nChan,nstart)=subplot('position',[X Y W H],'parent',Fig);
       
        
        
        %ExtractSignal
        if Channels(1,nChan)==Channels(2,nChan)
            BinData=ExtractContinuousData(Info.FilesDir,Info,Channels(1,nChan),RelStart,RelEnd);
            
            %conversion in Real unit
        
            VReal=ADC2Real(BinData,Info.Range(Channels(1,nChan))/2,Info.Gain(Channels(1,nChan)),Info.Offset(Channels(1,nChan)));
        else
            BinData1=ExtractContinuousData(Info.FilesDir,Info,Channels(1,nChan),RelStart,RelEnd);
            BinData2=ExtractContinuousData(Info.FilesDir,Info,Channels(2,nChan),RelStart,RelEnd);
            %conversion in Real unit
        
            VReal1=ADC2Real(BinData1,Info.Range(Channels(1,nChan))/2,Info.Gain(Channels(1,nChan)),Info.Offset(Channels(1,nChan)));
            VReal2=ADC2Real(BinData2,Info.Range(Channels(2,nChan))/2,Info.Gain(Channels(2,nChan)),Info.Offset(Channels(2,nChan)));
            VReal=VReal1-VReal2;
            
            clear VReal1 VReal2 BinData1 BinData2;
            
        end
        
        %add a DC remove
        if isnan(YZero(nChan))==1
            VReal=VReal-mean(VReal);
            OffsetY=0;
        else
            OffsetY=YZero(nChan);
        end 
               
        %filter
        
        if isempty(DHp)==0
           VReal=filtfilt(DHp.sosMatrix,DHp.ScaleValues,VReal);  
            
        end
        
        if isempty(DLp)==0
            VReal=filtfilt(DLp.sosMatrix,DLp.ScaleValues,VReal); 
        end
        
        
        Time=linspace(0,Duration,length(VReal));
        plot(h(nChan,nstart),Time,VReal,'color',CurrCol);
        
       
        set(h(nChan,nstart),'xlim',[0 Duration],'ylim',[-YScale(nChan)*2.5+OffsetY YScale(nChan)*2.5+OffsetY]);
        
        grid(h(nChan,nstart),'off');
        box(h(nChan,nstart),'off');
        axis(h(nChan,nstart),'off');
        
%         %add title and label
%         if nstart==1
%             ylabel(h(nChan,nstart),Info.ChLabel{nChan});
%         end
%         
%         
%         if nChan==nHaxe && isempty(LabelState)
%             title(h(nChan,nstart),LabelState{nstart})
%             
%             
%         end

        %add the line ref for scale
        if nstart==1 %&& nChan==1
            hold on;
            if nChan==1
            plot([0 XScale ],[-YScale(nChan)*2.5+OffsetY -YScale(nChan)*2.5+OffsetY],'k','linewidth',2);
            end
            plot([0 0],[-YScale(nChan)*2.5+OffsetY -YScale(nChan)*2.5+YScale(nChan)+OffsetY],'k','linewidth',2);
            
            
        end
         
        
        
    end
    
    linkaxes(h(nChan,:),'y');
    
end

% linkaxes(h,'x');

%create folder %s\\Computed Data
mkdir(sprintf('%s\\Computed Data',CurrDir));
mkdir(sprintf('%s\\Computed Data\\RawSignalFig\\',CurrDir));
FigFileName=sprintf('%s\\Computed Data\\RawSignalFig\\RawSignal_%s_%s.fig',CurrDir,Pref,Info.ExpFileName(1:end-4));
PdfFileName=sprintf('%s\\Computed Data\\RawSignalFig\\RawSignal_%s_%s.pdf',CurrDir,Pref,Info.ExpFileName(1:end-4));
MatFileName=sprintf('%s\\Computed Data\\RawSignalFig\\RawSignal_%s_%s.mat',CurrDir,Pref,Info.ExpFileName(1:end-4));

save(MatFileName);

hgsave(Fig, FigFileName, '-v7.3');
set(Fig,'PaperPositionMode','Auto','PaperOrientation','landscape');

print(Fig,PdfFileName,'-dpdf','-r600','-bestfit');

'done'
