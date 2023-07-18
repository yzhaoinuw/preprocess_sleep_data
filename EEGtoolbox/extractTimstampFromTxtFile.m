function [ImTimeStamp,ImNum,ImErr,meanfps,totalfps]=extractTimstampFromTxtFile(FileName,Type)

%if Type is 0 ImTimeStamp is in datenum format if Type is 1 ImTimeStamp is in datevec
if nargin==0
    FileName=[];
    Type=0; %datenum
elseif nargin==1
    Type=0; %datenum
end
    
if isempty(FileName)
    try
    load TempDir
    [f,DefaultDir]=uigetfile('*.txt','choose the Timer txt file',DefaultDir);
    catch
        [f,DefaultDir]=uigetfile('*.txt','choose the Timer txt file');
    end

end
    
 %import data from file   

%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: text (%q)
%	column4: text (%q)
%   column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%q%q%f%[^\n\r]';

%% Open the text file.
fileID = fopen(FileName,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
dataArray([1, 2, 5]) = cellfun(@(x) num2cell(x), dataArray([1, 2, 5]), 'UniformOutput', false);
Timercamera = [dataArray{1:end-1}];

% 
ImNum=[Timercamera{:,1}];
% ImTag=[Timercamera{:,2}];
% ImErr=[Timercamera{:,5}];

ImTimeStamp=Timercamera(:,3);
ImErr=[Timercamera{:,5}];
%remove err image
idx=ImErr==0;
idx(1:2)=0;

try
ImTimeStamp=datevec(ImTimeStamp,'yyyy-mm-dd_HH-MM-SS.FFF');
catch
ImTimeStamp=datevec(ImTimeStamp,'yyyy-mm-dd_HH-MM-SS-FFF');
end
%remove err image
idx=ImErr==0;
idx(1:2)=0;


%evaluate sliding fps
win=100; %windows in nb image

lim=1:100:sum(idx);

%Idok=ImNum(idx)+1;
Idok=1:length(ImTimeStamp);
for n=1:length(lim)-1
   
    id1=Idok(lim(n));
    id2=Idok(lim(n+1));
    
    
    try
    FpsIm(n)=(id2-id1+1)/etime(ImTimeStamp(id2,:),ImTimeStamp(id1,:));
    end
end



meanfps=mean(FpsIm);
totalfps=(Idok(end)-Idok(1)+1)/etime(ImTimeStamp(Idok(end),:),ImTimeStamp(Idok(1),:));
if Type==0
    ImTimeStamp=datenum(ImTimeStamp);
end





