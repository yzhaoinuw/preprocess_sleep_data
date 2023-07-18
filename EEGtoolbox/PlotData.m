function PlotData(MatFileName,Param,Statelist,StateColor)
%this function display the point cloud of data gerenated by extract parameters from slipanalysis

%matfile name coule be empty, if not this is the folder and filename of the file to open
%Param is a matrix containing two or three number that are the paremeter
%number in the DataProcess Matrix. ex: [2 4] or [2 3 8] 
%Statlist is a matrix that containt the statecode to plot ex  [1 2 4], could be empty
%StateColor is a matrix for each state to plot ex [0 0 1;1 0 0;0 1 0], coud be empty
%call Exemple PlotData('c:\Data\data.mat',[2 4 5]);   or   PlotData([],[1 4])  or
% or PlotData([],[1 4],[1 2 4],[0 0 1;1 0 0;0 1 0])   to the plot WK SWS PS with default color and code 

if isempty(MatFileName)==1
    uiload
end

Hypno=DataProcess(:,1);
Data=DataProcess(:,4:end);
if exist('Statelist')==0
    Statelist=unique(Hypno);
end

if exist('StateColor')==0
    StateColor=hsv(length(Statelist));
end

figure;

for n=1:length(Statelist)
    
    Id=Hypno==Statelist(n);
    
    if length(Param)==2
        plot(Data(Id,Param(1)),Data(Id,Param(2)),'color',StateColor(n,:),'marker','.','linestyle','none');hold on;
    elseif length(Param)==3
        plot3(Data(Id,Param(1)),Data(Id,Param(2)),Data(Id,Param(3)),'color',StateColor(n,:),'marker','.','linestyle','none');hold on;
    end
    
    grid on;
   
end
grid on;

xlabel(Labels(Param(1)+3));
ylabel(Labels(Param(2)+3));

xlim([quantile(Data(:,Param(1)),0.01) quantile(Data(:,Param(1)),0.99)]);
ylim([quantile(Data(:,Param(2)),0.01) quantile(Data(:,Param(2)),0.99)]);

if length(Param)==3
    zlabel(Labels(Param(3)+3));
    zlim([quantile(Data(:,Param(3)),0.01) quantile(Data(:,Param(3)),0.99)]);
end






