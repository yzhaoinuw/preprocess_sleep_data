function Acti=GetActimetry(Info,ActiName,T0,TEnd)

%This function get the actimetry from actiname between T0 and TEnd in min
%info come friom loadEXP and containt all the information about the file
if isempty(Info)==1
    Info=loadEXP;
end

AllActName={ Info.ActiFiles.Name};

CurrActi=find(strcmp(AllActName,ActiName));
DataActiTemp=[];
TimeActi=[];
DataActi=[];
if isempty(CurrActi)==0
    for nfile=1:length(Info.ActiFiles(CurrActi).Files)
        ActiFullFileName=fullfile(Info.ActiFiles(CurrActi).Files(nfile).Dir,Info.ActiFiles(CurrActi).Files(nfile).FileName);
        %Read the correct actimetry
        fidacti=fopen(ActiFullFileName,'r+');
        fseek(fidacti,0,'bof');
       
        try 
            DataActiTemp=fread(fidacti,[2,inf],'int32');
            TSTartActi=etime(datevec(Info.ActiFiles(CurrActi).Files(nfile).TStart),datevec(Info.BinFiles(1).TStart));
            %offsettime=etime(datevec(Info.ActiFiles(CurrActi).Files(nfile).TStart),datevec(Info.BinFiles(nhyp).TStart));

            TimeActi=[TimeActi DataActiTemp(1,:)./1000+TSTartActi];
            DataActi=[DataActi abs(DataActiTemp(2,:))];

        end
        fclose(fidacti);
    end
    
    %get the correct part
    IdTime=TimeActi>=T0*60 & TimeActi<=TEnd*60;
    Acti.Time=TimeActi(IdTime);
    Acti.Data=DataActi(IdTime);
    
%     figure;plot(Acti.TimeActi/60,Acti.DataActi,'b');
%     xlabel('time (min)');
%     ylabel(ActiName);
%     grid on;
%     
    
    
else
    sprintf('Actimetry %s not found',ActiName);
end

