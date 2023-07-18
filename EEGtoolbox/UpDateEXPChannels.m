function s=UpDateEXPChannels(s,Id,Name,Offset,Gain,AcquisitionRangeMax,RecordType)

%s= xml2struct(XMLfilename);
if isfield(s.Animal.Acquisition.Channels,'Channel')==0
    NbChan=1;
else
    NbChan=length(s.Animal.Acquisition.Channels.Channel)+1;
end

s.Animal.Acquisition.Channels.Channel{NbChan}.Id=num2str(Id);
s.Animal.Acquisition.Channels.Channel{NbChan}.Enable='true';
s.Animal.Acquisition.Channels.Channel{NbChan}.Name=Name;
s.Animal.Acquisition.Channels.Channel{NbChan}.Offset=num2str(Offset);
s.Animal.Acquisition.Channels.Channel{NbChan}.Gain=num2str(Gain);
s.Animal.Acquisition.Channels.Channel{NbChan}.AcquisitionRangeMax=num2str(AcquisitionRangeMax);
s.Animal.Acquisition.Channels.Channel{NbChan}.RecordType=num2str(RecordType);

s.Animal.Acquisition.NbChan=num2str(NbChan);

%struct2xml(s,XMLfilename);