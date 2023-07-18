function Info=renameAct(d)
%change the time act file with the closest time avi
try 
    load tempDir;
catch
    DefaultDir='c:\';
end;

if DefaultDir==0
    DefaultDir='c:\';
end

if nargin==0
    [d]=uigetdir(DefaultDir);
    DefaultDir=d;
    save TempDir DefaultDir;
end

filelist=dir(d);
nact=0;
navi=0;
actifile={};
avifile={};
for n=1:length(filelist)
    currFile=filelist(n).name;
    if length(currFile)>5 && strcmp(currFile(end-2:end),'act')==1
        nact=nact+1;
        actifile{nact}=currFile;
        idlim=find(currFile=='_',1,'last');
        actiname{nact}=currFile(idlim+1:end-4);
        actitime{nact}=datenum(currFile(idlim-23:idlim-1),'yyyy-mm-dd_HH-MM-SS-FFF');
        
    elseif length(currFile)>5 && strcmp(currFile(end-2:end),'avi')==1
        navi=navi+1;
        idlim=find(currFile=='_',1,'last');
        avifile{navi}=currFile;
        aviname{navi}=currFile(idlim+1:end-4);
        avitime{navi}=datenum(currFile(idlim-23:idlim-1),'yyyy-mm-dd_HH-MM-SS-FFF');
    end
    
end


%rename actifile with the correct avi file 
Info={};
for nact=1:length(actifile)

    idavi=find((strcmp(aviname,actiname{nact}) &  abs(actitime{nact}-cell2mat(avitime)).*3600*24<0.5 & abs(actitime{nact}-cell2mat(avitime)).*3600*24>=0)==1);
    Info{nact,1}=actifile{nact};
    if length(idavi)==1
        
        difftime=abs(actitime{nact}-avitime{idavi}).*3600*24;

        Info{nact,2}=avifile{idavi};
        Info{nact,4}=difftime;
        
        if difftime==0
            Info{nact,5}='Time already Ok';    
            Info{nact,3}='file not renamed'
        else
            Info{nact,5}='Ok';
            newname=[avifile{idavi}(1:end-3) 'act'];
            Info{nact,3}=newname;

           [status,message]=copyfile(fullfile(d,actifile{nact}),fullfile(d,newname));
           delete(fullfile(d,actifile{nact}));
        end
        
    elseif length(idavi)==0
        Info{nact,5}='no avifile found';
    else
        Info{nact,5}='multiple avi file detected';
    end
        
end
    
    
    
'Done'

