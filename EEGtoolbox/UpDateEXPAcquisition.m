function s=UpDateEXPAcquisition(s,NbChan,SamplingRate)

%s= xml2struct(XMLfilename);

s.Animal.Acquisition.NbChan=num2str(NbChan);
s.Animal.Acquisition.SamplingRate=num2str(SamplingRate,20);


%struct2xml(s,XMLfilename);