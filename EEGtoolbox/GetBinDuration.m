function NbSample=GetBinDuration(filename,nbchan)
fid=fopen(filename);
if fid~=-1
    fseek(fid,0,'eof');
    NbSample=ftell(fid)/2/nbchan;
    fclose(fid);
else
    NbSample=0;
end