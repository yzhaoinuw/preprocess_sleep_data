function [meanDataPerSeq,TimeLim]=MeanDataperseq(DataFileName,Pref,SeqDuration,State)
%this function compute the mean of the parameters (meanDataPerSeq) evaluated from
%slipanalysis (compute feature) per sequence of SeqDuration (s)
%DataFileName is the mat file name (with directory)  that contain the variable DataProcess came from compute feature.m could be empty
%Pref prefix for the output files
% State is a matrix that could be empty to filter the states
%Timelim contain the time line in min for the computed data (start like the
%time data in dataprocess from the first data aquired in the first binary
%file of the exp file.

if isempty(DataFileName) ==1
	[f,d]=uigetfile('*.mat','pick the mat file with the data');
	DataFileName=fullfile(d,f);
	
end
load(DataFileName);
if exist('DataProcess')==1

	[directory,Filename,ext]=fileparts(DataFileName);
      
    FigOutfilename=fullfile(directory,[Pref '-' Filename '.fig']);
    XlsOutfilename=fullfile(directory,[Pref '-' Filename '.xlsx']);
    MatOutfilename=fullfile(directory,[Pref '-' Filename '.mat']);
    
    
	Time=DataProcess(:,2);
	nbseq=floor((Time(end)-Time(1))/SeqDuration);
	nbdata=find((Time-Time(1))<SeqDuration,1,'last');

	%filter state
	idstate=ismember(DataProcess(:,1),State);
	DataProcess(~idstate,4:end)=NaN;
	idmax=nbseq*nbdata;
    
    SeqDuration=SeqDuration/60;
	TimeLim=Time(1)/60+linspace(SeqDuration/2,SeqDuration*nbseq-SeqDuration/2,nbseq);
    
    Labels=Labels(4:end);
    Out=[{'Time (Hours)'} num2cell(TimeLim)];
        
    meanDataPerSeq=[];
	for np=4:size(DataProcess,2)
		TabSeq=reshape(DataProcess(1:idmax,np),nbdata,nbseq);
		meanDataPerSeq(np-3,:)=nanmean(TabSeq,1);
        Out{np-2,1}=Labels{np-3};
        Out(np-2,2:1+size(meanDataPerSeq,2))=num2cell(meanDataPerSeq(np-3,:));
	end
	
    
    fig=figure;
    col=hsv(size(DataProcess,2)-3);
    for n=1:size(DataProcess,2)-3
        plot(TimeLim,meanDataPerSeq(n,:),'color',col(n,:));
        hold on;       
    end
    grid on;
    xlabel('time min');
   
    legend(Labels);
    
    
    
    
    save(MatOutfilename,'TimeLim','meanDataPerSeq','Labels');
    xlswrite(XlsOutfilename,Out);
    hgsave(fig, FigOutfilename, '-v7.3')
	
else
	errordlg('the file does not contain the matrix DataProcess');
end	


