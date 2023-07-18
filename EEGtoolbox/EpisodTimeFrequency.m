function [TF,F]=EpisodTimeFrequency(EpisodLim,Info,ChanNum,MTT)

TF=[];
MTT.params.Fs=Info.Fs;
for nepisod=1:length(EpisodLim.Deb)


    %exctract the uint16 raw data from bin file


    if size(ChanNum,2)==1% if no substraction
        if size(ChanNum,1)==1
            ChanNum(2,1)=ChanNum(1,1);
        end
    end

    if ChanNum(1,1)==ChanNum(2,1) %if no substraction
          BinData=ExtractContinuousData(Info.FilesDir,Info,ChanNum(1,1),EpisodLim.Deb(nepisod),EpisodLim.End(nepisod));

         %conversion, in real unit
         DataChan1=ADC2Real(BinData,Info.Range(ChanNum(1,1))/2,Info.Gain(ChanNum(1,1)),Info.Offset(ChanNum(1,1)));
        clear BinData;
      %  ChanName{1}= Info.ChLabel{ChanNum(1,1)};

    elseif ChanNum(1,1)~=ChanNum(2,1) %if substraction 


        BinData1=ExtractContinuousData(Info.FilesDir,Info,ChanNum(1,1),EpisodLim.Deb(nepisod),EpisodLim.End(nepisod));

        BinData2=ExtractContinuousData(Info.FilesDir,Info,ChanNum(2,1),EpisodLim.Deb(nepisod),EpisodLim.End(nepisod));

        %conversion, in real unit
        VReal1=ADC2Real(BinData1,Info.Range(ChanNum(1,1))/2,Info.Gain(ChanNum(1,1)),Info.Offset(ChanNum(1,1)));

        %conversion, in real unit
        VReal2=ADC2Real(BinData2,Info.Range(ChanNum(2,1))/2,Info.Gain(ChanNum(2,1)),Info.Offset(ChanNum(2,1)));

        clear BinData1 BinData2;

        DataChan1=VReal1-VReal2;

        clear  VReal1 VReal2;

       % ChanName{1}= [Info.ChLabel{ChanNum(1,1)} '-' Info.ChLabel{ChanNum(2,1)}];

    end
    
    [S,~,F]=mtspecgramc(DataChan1,MTT.Movingwin,MTT.params);

    TF=[TF;S];
end
