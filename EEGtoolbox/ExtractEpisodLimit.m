function StateLim=ExtractEpisodLimit(Hypno)



allstate=unique(Hypno); 
Time=0:length(Hypno)-1;

for nState=1:length(allstate)
    Code=allstate(nState);

    id=Hypno==Code;
    StateLim(nState).Code=Code;
    StateLim(nState).Deb=Time(diff([0 id'])==1);
    StateLim(nState).End=Time(diff([id' 0])==-1)+1;
    StateLim(nState).Duration=StateLim(nState).End-StateLim(nState).Deb;

end

 

 
 
 