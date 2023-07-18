function VReal=ADC2Real(ADC,MaxRange,Gain,Offset)
%MaxRange half whole range
VReal=(ADC.*repmat(MaxRange',1,size(ADC,2)).*2./2^16-repmat(Offset',1,size(ADC,2)))./repmat(Gain',1,size(ADC,2));


