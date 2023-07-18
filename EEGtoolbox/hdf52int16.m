function hdf52int16(filename,Gain,Range)

if isempty(filename)==1
    [f,d]=uigetfile('*.h5');
    filename=fullfile(d,f);
end

fid = H5F.open(filename);
dset_id = H5D.open(fid,'/raw/samples');
data = H5D.read(dset_id);
H5F.close(fid);


%conversion in binary file
MaxRange=Range/2;
Offset=0;
type='int16';
ADC=Real2ADC(data,MaxRange,Gain,Offset,type);

[d,f,ext]=fileparts(filename)
binfilename=fullfile(d,[f '.bin']);

fidbin=fopen(binfilename,'w+');
fwrite(fidbin,ADC,type);
fclose(fidbin);

'done'
