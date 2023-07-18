function BinData=removesat(BinData,Precision,SatPer)

%this function remove saturation by an iunterpolation with the nearest neighbor
%bindata ais the data to proceed
%Precision is uint16 or int16
%SatPer is the % of the range that is consider as a saturation

%remove sat (oneiros)
if strcmpi(Precision,'uint16')
    minADC=2^15-SatPer*2^16/2/100;
    maxADC=2^15+SatPer*2^16/2/100;
else
    minADC=-SatPer*2^16/2/100;
    maxADC=SatPer*2^16/2/100; 
end



for n=1:size(BinData,1)
    idsat=BinData(n,:)>=maxADC | BinData(n,:)<=minADC;    

    if sum(idsat)>0
        X1=1:size(idsat,2);
        X0=X1(~idsat(:));
        Y0=BinData(n,:);
        Y0(idsat)=[];

        BinData(n,:) = interp1(X0,Y0,X1,'nearest');
    end
  
end 