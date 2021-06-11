function h = MyScatterBarPlotSD(y,x,barWidth,xColor)
%MyScatterBarPlot Plot scatter plot and Bar plot with SD error bar
%
% function MyScatterBarPlot(data,group,groupColor)
%
%
% Inputs
% y - value
% x - group variable for y
%
% Dohoung Kim - June 2015
group = unique(x);
nGroup = length(group);

hold on;
for iGroup = 1:nGroup
    yPoint = y(x==group(iGroup) & ~isnan(y));
    nPoint = sum(~isnan(yPoint));
    xPoint = iGroup + 0.75*barWidth*(rand(nPoint,1)-0.5);
       
    yMean = nansum(yPoint)/nPoint;
    ySem = nanstd(yPoint);
    
    h(iGroup).bar = bar(iGroup,yMean,'FaceColor',xColor{iGroup},'LineStyle','none','BarWidth',barWidth);
    h(iGroup).errorbar = errorbar(iGroup,yMean,ySem,'LineWidth',2,'Color',xColor{iGroup});
    errorbarT(h(iGroup).errorbar,0.2,0.5);

        
    plot(xPoint,yPoint,'LineStyle','none','Marker','.','MarkerSize',6,'Color',xColor{iGroup});
    
end
set(gca,'Box','off','TickDir','out','FontSize',4,'LineWidth',0.35,...
    'XLim',[0.5 nGroup+0.5],'XTick',1:nGroup);