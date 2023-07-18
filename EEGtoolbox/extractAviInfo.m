function Info=extractAviInfo(filename)

fid=fopen(filename);
fseek(fid,0,'bof');


%search fps
STRH=1752331379;
maxpos=30*4;
OffsetSTRH=6*4;
RawV=searchV(OffsetSTRH,fid,STRH,4);
Info.Fps=RawV(2)/RawV(1);
Info.Start=RawV(3);
Info.NumFrame=RawV(4);


%search NumFrame

AVIH=1751742049;
maxpos=30*4;
OffsetAVIH=5*4;
Raw=searchV(OffsetAVIH,fid,AVIH,6);

Info.Width=Raw(5);
Info.Height=Raw(6);

fclose(fid);


%search Resolution

function Vout=searchV(Offset,fid,tag,size)


maxpos=30*4;
pos=0;
fseek(fid,0,'bof');
CurrV=0;

while CurrV~=tag && pos<=maxpos
    CurrV=(fread(fid,1,'uint'));
    pos=ftell(fid);
end

if pos==maxpos %not detected
    Vout=[];   
else
    fseek(fid,Offset,'cof');
    Vout=fread(fid,size,'uint');
end

