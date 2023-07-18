function ADC=Real2ADC(VReal,HalfRange,Gain,Offset,type)
if nargin<=4
    type=uint16;
end
if strcmpi(type,'int16')
    ADC=int16((VReal.*repmat(Gain',1,size(VReal,2))+repmat(Offset',1,size(VReal,2))).*2^16./repmat(HalfRange',1,size(VReal,2))./2);
elseif strcmpi(type,'uint16')
    ADC=uint16((VReal.*repmat(Gain',1,size(VReal,2))+repmat(Offset',1,size(VReal,2))).*2^16./repmat(HalfRange',1,size(VReal,2))./2);
end
