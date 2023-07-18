function ChangeAviFrameRate(FileName,NewFrameRate)

NumFR= uint32(NewFrameRate*10000000)
DeNumFR=  uint32(10000000)
fid=fopen(FileName,'r+');
fseek(fid,0,'bof');

try
    %search fps
    STRH=1752331379;
    maxpos=30*4;
    OffsetSTRH=6*4;

    size=2;
    pos=0;
    CurrV=0;

    while CurrV~=STRH && pos<=maxpos
        CurrV=(fread(fid,1,'uint'));
        pos=ftell(fid);
    end

    if pos==maxpos %not detected
        Vout=[];   
    else
        fseek(fid,OffsetSTRH,'cof');
        fwrite(fid,DeNumFR,'uint');
        fwrite(fid,NumFR,'uint');
        fseek(fid,-2*4,'cof');
        RawV=fread(fid,size,'uint');

        FpsInFile=RawV(2)/RawV(1)
    end
end
fclose(fid);