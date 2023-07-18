function xmlFileName=CreateEXPfile(filename,AnimalName,FileOrigin,StimulationChannelId,Precision,StateList,EpochDur)
if nargin<=6
    EpochDur=5;
end
docNode = com.mathworks.xml.XMLUtils.createDocument('Animal');
docRootNode = docNode.getDocumentElement;
%docRootNode.setAttribute('attr_name','attr_value');
%init animal field
    thisEnabled = docNode.createElement('Enabled'); 
    thisEnabled.appendChild... 
        (docNode.createTextNode('true'));
    docRootNode.appendChild(thisEnabled);
    
    thisName = docNode.createElement('Name'); 
    thisName.appendChild... 
        (docNode.createTextNode(AnimalName));
    docRootNode.appendChild(thisName);
    
    thisFileOrigin = docNode.createElement('FileOrigin'); 
    thisFileOrigin.appendChild... 
        (docNode.createTextNode(FileOrigin));
    docRootNode.appendChild(thisFileOrigin);
    
    thisAcquisition = docNode.createElement('Acquisition'); 
    
        thisSamplingRate = docNode.createElement('SamplingRate');      
         thisAcquisition.appendChild(thisSamplingRate);
         
         thisEnable = docNode.createElement('Enable');      
         thisEnable.appendChild... 
        (docNode.createTextNode('true'));
         thisAcquisition.appendChild(thisEnable);
         
         thisFiles = docNode.createElement('Files');      
         thisAcquisition.appendChild(thisFiles);
         
         thisChannels = docNode.createElement('Channels');      
         thisAcquisition.appendChild(thisChannels);
         
         thisFiles = docNode.createElement('NbChan');      
         thisAcquisition.appendChild(thisFiles);
         
         thisHeaderOffset = docNode.createElement('HeaderOffset');   
         thisHeaderOffset.appendChild... 
        (docNode.createTextNode('0'));
         thisAcquisition.appendChild(thisHeaderOffset);
         
         thisAcquisitionType = docNode.createElement('AcquisitionType');
         thisAcquisitionType.appendChild... 
        (docNode.createTextNode(Precision));
         thisAcquisition.appendChild(thisAcquisitionType);
        
        
    docRootNode.appendChild(thisAcquisition);
    
    thisElement = docNode.createElement('Videos'); 
    
    docRootNode.appendChild(thisElement);
    
    thisHypno = docNode.createElement('Hypnogram'); 
         thisHypSampling = docNode.createElement('SamplingRate');   
         thisHypSampling.appendChild... 
        (docNode.createTextNode(num2str(1/EpochDur)));
         thisHypno.appendChild(thisHypSampling);
         
         thisHypEnable = docNode.createElement('Enable');   
         thisHypEnable.appendChild... 
        (docNode.createTextNode('true'));
         thisHypno.appendChild(thisHypEnable);
         
         thisHypFiles = docNode.createElement('Files');   
         thisHypno.appendChild(thisHypFiles);
         
          thisEEGChannel = docNode.createElement('EEGChannel');   
         thisHypno.appendChild(thisEEGChannel);
         
          thisEMGChannel = docNode.createElement('EMGChannel');   
         thisHypno.appendChild(thisEMGChannel);
         
          thisConfigFile = docNode.createElement('ConfigFile');   
         thisHypno.appendChild(thisConfigFile);
         
          thisStatuses = docNode.createElement('Statuses');
         
          
          if nargin==5
                
              Level=[0 1 2 3 4 5];
              Label={'ND' 'WK' 'SWS' 'TR' 'PS' 'Artf'};
              Key={0 1 2 3 4 5};
              Color=[4286611584 4278190335 4294901760 4294967040 4278222848 4278190080];
          elseif nargin>=6
              
              CodeList=[StateList(:).Code];
              [~,id]=sort(CodeList);
              
              StateList=StateList(id);
              
              
              for nState=1:length(StateList)
                  Level(nState)=StateList(nState).Code;
                  Label{nState}=StateList(nState).Label;
                  CurrKey=StateList(nState).Key;
                  if isempty(strfind(StateList(nState).Key,'numpad'))==0
                      Key{nState}=CurrKey(7:end);
                  else
                      Key{nState}=CurrKey;
                  end
                  Color(nState)=RGB2Dec(StateList(nState).Color);
              end
              
              
          end 
          
          for nstatus=1:length(Level)
             thisS = docNode.createElement('Status');
                
                thisS2 = docNode.createElement('Level');
                thisS2.appendChild(docNode.createTextNode(num2str(Level(nstatus))));
                thisS.appendChild(thisS2);
                
                thisS2 = docNode.createElement('Label');
                thisS2.appendChild(docNode.createTextNode(Label(nstatus)));
                thisS.appendChild(thisS2);
                
                thisS2 = docNode.createElement('Key');
                thisS2.appendChild(docNode.createTextNode(Key(nstatus)));
                thisS.appendChild(thisS2);
                
                thisS2 = docNode.createElement('Color');
                thisS2.appendChild(docNode.createTextNode(num2str(Color(nstatus))));
                thisS.appendChild(thisS2);
                
                
             thisStatuses.appendChild(thisS);
          end
        
       thisHypno.appendChild(thisStatuses);
         
   
    docRootNode.appendChild(thisHypno);

    thisElement = docNode.createElement('SleepDeprivation'); 
   
    docRootNode.appendChild(thisElement);
    
    thisElement = docNode.createElement('StimulationChannelId'); 
    thisElement.appendChild... 
        (docNode.createTextNode(StimulationChannelId));
    docRootNode.appendChild(thisElement);
    

xmlFileName = fullfile(FileOrigin,[filename '.xml']);
xmlwrite(xmlFileName,docNode);
