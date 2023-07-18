function CallMeanPowerSpectrumPerGroupExemple
%PowerSpectrum Setting
Settings.MTT.Movingwin=[3 0.5];%win and step
Settings.MTT.BandWidth=2;
Settings.MTT.Taper=5;

Settings.MTT.params.trialave = 1;
Settings.MTT.params.Fs = [];%will be filled after
Settings.MTT.params.tapers = [Settings.MTT.BandWidth Settings.MTT.Movingwin(1) 2*Settings.MTT.BandWidth*Settings.MTT.Movingwin(1)-Settings.MTT.Taper];
Settings.MTT.params.fpass = [0.5 45];
Settings.MTT.params.pad = 1;

Settings.Filt=[];
Settings.ApplyICA=0;
Settings.UseWhitening=0;
Settings.ARmodel=[];
Settings.EpochDuration=5;
Settings.FigSave=0;
Settings.LightOn=[];
Settings.ProcessDataPerXMin=[];
Settings.Artef=[];

Settings.Band2Analyse={[1 4] [4 9]};

Analyse(1).Mice={'S1' 'S2' 'S3' 'S4' 'S6' 'S7' 'S8' 'S9' 'S10' 'S11' 'S12'};
Analyse(1).Baseline='BL (10H-16H)-SWS';
Analyse(1).Manip='APSD (10H-16H)-SWS';
Analyse(1).Duration=360;
Analyse(1).TransitionTime2Remove=10;

Analyse(2).Mice={'S1' 'S2' 'S3' 'S4' 'S6' 'S7' 'S8' 'S9' 'S10' 'S11' 'S12'};
Analyse(2).Baseline='BL (16H-20H)-SWS';
Analyse(2).Manip='PSR (16H-20H)-SWS';
Analyse(2).Duration=240;
Analyse(2).TransitionTime2Remove=10;

Analyse(3).Mice={'S1' 'S2' 'S3' 'S4' 'S6' 'S7' 'S8' 'S9' 'S10' 'S11' 'S12'};
Analyse(3).Baseline='BL (16H-20H)-PS';
Analyse(3).Manip='PSR (16H-20H)-PS';
Analyse(3).Duration=240;
Analyse(3).TransitionTime2Remove=10;

%Settings extraction states
Settings.MinEpisodDuration=0;%0s  %5s in min = 0.083333333333333
Settings.ArtefactDetection=0;

SaveTF=0;

MeanPowerSpectrumPerGroup(Analyse,Settings,SaveTF)