
OutNamePref='test'; %nom du fichier de sortie
T0=0;% start time in min
TotDuration=20;% durée d'analyse en min

Info=loadEXP; %load the exp file to process

Pass=1; %for Time frequency process

Artef=[]; %no artefact used

FigSave=1; %save fig in tif and fig format

%define the epoch size and the number of epoch transition to remove
BoutSize=5;
NBTranstionWin2Remove=0; %no transition to remove

ProcessDataPerXMin=[];%the analyse will not be cut

%filter design

    %generate filters
    %filter band pass ([10 40 hz]) order 10
    Filt(1).hp  = fdesign.bandpass('N,Fc1,Fc2',10,10,40,Info.Fs);
    Filt(1).Hd = design(Filt(1).hp,'butter');
    
    %high pass filter at 0.2hz order 10
    Filt(2).hp  = fdesign.highpass('N,Fc', 10,0.2,Info.Fs);
    Filt(2).Hd = design(Filt(2).hp,'butter');

    %high pass filter at 20hz order 2
    Filt(3).hp  = fdesign.highpass('N,Fc', 2,20,Info.Fs);
    Filt(3).Hd = design(Filt(3).hp,'butter');
    
    %low pass filter at 5hz order 2
    Filt(4).hp  = fdesign.lowpass('N,Fc', 2,5,Info.Fs);
    Filt(4).Hd = design(Filt(4).hp,'butter');
    
    

%multitaper design

    %generate parameter for Multitaper
    MTT(1).Movingwin=[2 0.5];%win and step
    MTT(1).BandWidth=2;
    MTT(1).Taper=3;

    MTT(1).params.trialave = 1;
    MTT(1).params.Fs = [];%will be filled after
    MTT(1).params.tapers = [MTT(1).BandWidth MTT(1).Movingwin(1) 2*MTT(1).BandWidth*MTT(1).Movingwin(1)-MTT(1).Taper];
    MTT(1).params.fpass = [0 20];
    MTT(1).params.pad = 1;
    
    MTT(2).Movingwin=[5 1];%win and step
    MTT(2).BandWidth=2;
    MTT(2).Taper=10;

    MTT(2).params.trialave = 1;
    MTT(2).params.Fs = [];%will be filled after
    MTT(2).params.tapers = [MTT(2).BandWidth MTT(2).Movingwin(1) 2*MTT(2).BandWidth*MTT(2).Movingwin(1)-MTT(2).Taper];
    MTT(2).params.fpass = [0 40];
    MTT(2).params.pad = 1;
    
    
    %very low frequency
    MTT(3).Movingwin=[2 0.1];%win and step
    MTT(3).BandWidth=25;
    MTT(3).Taper=5;

    MTT(3).params.trialave = 1;
    MTT(3).params.Fs = [];%will be filled after
    MTT(3).params.tapers = [MTT(3).BandWidth MTT(3).Movingwin(1) 2*MTT(3).BandWidth*MTT(3).Movingwin(1)-MTT(3).Taper];
    MTT(3).params.fpass = [0 30];
    MTT(3).params.pad = 1;
    
%define the process

    %process is define a describe here
    %first cell is the channel if 2 value in the vector (ex: '[2 3]') means a difference
    %second cell is the kind of parameter to evaluate
    %   1 Compute TimeFrequency with MTT
    %   2 is heart rate  
    %   3 raw signal
    %   4 Peak frequency in specfic band
    %   5 Band power
    %   6 pattern correlation
    %   7 Band ratio
    %   8 Peak detection


    %the third cell is the option for each parameter
    %   1: MTT parameters number
    %   2: [minpeakdist (in s) minpeakvalue]
    %   3: [] 
    %   4: [fmin1 fmax1]; use only if Time frequency have been calculated
    %   5: [fmin1 fmax1]; use only if Time frequency have been calculated
    %   6: mat file which contain the reference pattern in a matrix named
    %   Pattern with a [m,1] size the mat file should be in the same  directory of the exp file pattern is in uint16 value
    %   7: [fmin1 fmax1;fmin2 fmax2] use only if Time frequency have been calculated
    %   8: [MinPeakDist MinPeakThreshold] MinPeakDist is the minimum time in s between two peaks,MinPeakThreshold if thresh is negative find negative peak, MinPeakThreshold is in number of std from the mean

    %the fourth cell is the decimate factor %if empty mean no subsampling

    %the fifth cell is the number of filter to use, if empty no filter else
    %should be 1 or 2 depending on filter designed above or [] if no filter

    %the sixth cell is the ylim for display or clim depending on the kind
    %of analysis, if empty limit is set automatically
    
    %the seventh cell is if yes or no artefact removing is used  (base on ICA detection);



ProcessInfo(1,1:7)={[8 3] 1 1 2 [] [-211.551267250141 -78.5791906280928] 0}; %calcul a realiser
ProcessInfo(2,1:7)={[9 3] 1 1 2 [] [] 0}; %calcul a realiser
ProcessInfo(3,1:7)={[10 3] 1 1 2 [] [] 0}; %calcul a realiser
ProcessInfo(4,1:7)={[11 3] 1 1 2 [] [] 0}; %calcul a realiser
ProcessInfo(5,1:7)={[14 3] 1 1 2 [] [] 0}; %calcul a realiser
ProcessInfo(6,1:7)={[15 3] 1 1 2 [] [] 0}; %calcul a realiser
ProcessInfo(7,1:7)={[16 3] 1 1 2 [] [] 0}; %calcul a realiser

    
    
%call the computeWUW script

ComputeWUW(OutNamePref,T0,TotDuration,ProcessInfo,Pass,Info,Artef,FigSave,...
    BoutSize,NBTranstionWin2Remove,ProcessDataPerXMin,Filt,MTT)