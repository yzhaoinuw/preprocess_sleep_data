function Bin2EDF(ExpFileName,Mutiplefile)
%this function convert EXP file with binary file into edf file 
%ExpFileName is the full name including the dircetory and filename of the exp file
%MutlipleFile is an option to generate one file per binary file if value is 1
try
    if nargin==0 ||  isempty(ExpFileName)
        [f,d]=uigetfile('*.exp');
        ExpFileName=fullfile(d,f);  
    end


    [Directory,Prefname]=fileparts(ExpFileName);

    if nargin==0 || nargin==1 || isempty(Mutiplefile) 
        Mutiplefile=1;
    end
    barre=waitbar(0,'conversion in progress');
    Info=loadEXP(ExpFileName,'No');

    NbBinfiles=length(Info.BinFiles);
    FirstStart=Info.BinFiles(1).TStart;

    Header = struct();
    Header.samplingrate = Info.Fs;
    Header.channels     =  Info.ChLabel;


    Header.subject.ID    = 'X';
    Header.subject.sex   = 'X';
    Header.subject.name  = Info.Name;
    Header.subject.year  = NaN;
    Header.subject.month = NaN;
    Header.subject.day   = NaN;

    
    if Mutiplefile==1 || NbBinfiles==1

         for nfile=1:NbBinfiles
             waitbar(nfile/NbBinfiles,barre);
            CurrentStart=Info.BinFiles(nfile).TStart;
            EdfFilename=[Prefname '_' datestr(CurrentStart,'yyyy-mm-dd_HH-MM-SS-FFF')];

            CurrVecDate= datevec(CurrentStart);
            Header.year         = CurrVecDate(1);
            Header.month        = CurrVecDate(2);
            Header.day          = CurrVecDate(3);
            Header.hour         = CurrVecDate(4);
            Header.minute       = CurrVecDate(5);
            Header.second       = CurrVecDate(6);
            Data=[];
            fid=fopen(fullfile(Info.BinFiles(nfile).Dir,Info.BinFiles(nfile).FileName),'r+');
            Data=ReadBinFile(Info.Fs,0,inf,fid,Info.NbRecChan,0,Info.Precision);
            fclose(fid);
            Data=ADC2Real(Data,Info.Range/2,Info.Gain,Info.Offset);
            lab_write_edf(fullfile(Info.FilesDir,EdfFilename), Data, Header);


        end




    else
        EdfFilename=[Prefname '_' datestr(FirstStart,'yyyy-mm-dd_HH-MM-SS-FFF')];
        CurrVecDate= datevec(FirstStart);
        Header.year         = CurrVecDate(1);
        Header.month        = CurrVecDate(2);
        Header.day          = CurrVecDate(3);
        Header.hour         = CurrVecDate(4);
        Header.minute       = CurrVecDate(5);
        Header.second       = CurrVecDate(6);

        Data=ExtractContinuousData(Info.FilesDir,Info,[],0,inf,[]);
        Data=ADC2Real(Data,Info.Range/2,Info.Gain,Info.Offset);
        lab_write_edf(fullfile(Info.FilesDir,EdfFilename), Data, Header);
    end

    warndlg('Conversion done')
catch
    errordlg('Conversion failed')
end
close barre;