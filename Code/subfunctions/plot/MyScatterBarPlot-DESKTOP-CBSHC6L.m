function h = MyScatterBarPlot(y,x,barWidth,yLim,xColor,sigGroup)
%MyScatterBarPlot Plot scatter plot and Bar plot with sem error bar
%
% function MyScatterBarPlot(data,group,groupColor)
%
%
% Inputs
% y - value
% x - group variable for y
%
%   Only compatible with <2014a version
%
% Dohoung Kim - June 2015
% Revision: June 2016
%   - Add significant line
group = unique(x);
nGroup = length(group);

hold on;
yPeak = zeros(1, nGroup);
yMax = zeros(1, nGroup);
yMin = zeros(1, nGroup);
for iGroup = 1:nGroup
    yPoint = y(x==group(iGroup) & ~isnan(y));
    nPoint = sum(~isnan(yPoint));
    xPoint = iGroup + 0.75*barWidth*(rand(nPoint,1)-0.5);
       
    yMean = nansum(yPoint)/nPoint;
    ySem = nanstd(yPoint)/sqrt(nPoint);
    
    yPeak(iGroup) = yMean+ySem;
    yMax(iGroup) = max([max(yPoint) yMean+ySem]);
    yMin(iGroup) = min([min(yPoint) yMean-ySem]);

    h(iGroup).bar = bar(iGroup,yMean,'FaceColor',xColor{iGroup},'LineStyle','none','BarWidth',barWidth);
    h(iGroup).errorbar = errorbar(iGroup,yMean,ySem,'LineWidth',2,'Color',xColor{iGroup});
    errorbarT(h(iGroup).errorbar, 0.2, 0.5);
    
    yPoint(yPoint>yLim(end)) = yLim(end);
    yPoint(yPoint<yLim(1)) = yLim(1);
    plot(xPoint,yPoint,'LineStyle','none','Marker','o','MarkerSize',2,...
        'MarkerFaceColor',xColor{iGroup},'MarkerEdgeColor','none');
end
yMax = max(yMax);
yMin = min(yMin);
if min(yMin)>=0; yMin=0; end;

if nargin==6 && iscell(sigGroup)
    nSig = length(sigGroup);
    for iS = 1:nSig
        if length(sigGroup{iS})==2
            if sum(isinteger(sigGroup{iS}))==2
            yHor = max(yPeak(sigGroup{iS}))+(yMax-yMin)*0.05;
            
            plot([sigGroup{iS}(1) sigGroup{iS}(1)], [yPeak(sigGroup{iS}(1))+(yMax-yMin)*0.025 yHor], ...
                'LineWidth', 0.35, 'Color', 'k');
            plot([sigGroup{iS}(2) sigGroup{iS}(2)], [yPeak(sigGroup{iS}(2))+(yMax-yMin)*0.025 yHor], ...
                'LineWidth', 0.35, 'Color', 'k');
            plot([sigGroup{iS}(1) sigGroup{iS}(2)], [yHor yHor], ...
                'LineWidth', 0.35, 'Color', 'k');
            text(mean(sigGroup{iS}), yHor+(yMax-yMin)*0.025, '*', ...
                'FontSize', 10, 'HorizontalAlignment', 'center');
            else
                yHor = max(yPeak([ceil(sigGroup{iS}),floor(sigGroup{iS})]))+(yMax-yMin)*0.05;
            
            plot([sigGroup{iS}(1) sigGroup{iS}(1)], [yPeak(nanmean(floor(sigGroup{iS}(1)),ceil(sigGroup{iS}(1))))+...
                (yMax-yMin)*0.025 yHor], ...
                'LineWidth', 0.35, 'Color', 'k');
            plot([sigGroup{iS}(2) sigGroup{iS}(2)], [yPeak(nanmean(floor(sigGroup{iS}(2)),ceil(sigGroup{iS}(2))))+...
                (yMax-yMin)*0.025 yHor], ...
                'LineWidth', 0.35, 'Color', 'k');
            plot([sigGroup{iS}(1) sigGroup{iS}(2)], [yHor yHor], ...
                'LineWidth', 0.35, 'Color', 'k');
            text(mean(sigGroup{iS}), yHor+(yMax-yMin)*0.025, '*', ...
                'FontSize', 10, 'HorizontalAlignment', 'center');
            end
        end
    end
end
    
set(gca,'Box','off','TickDir','out','FontSize',4,'LineWidth',0.35,...
    'XLim',[0.5 nGroup+0.5],'XTick',1:nGroup, ...
    'YLim', [yMin yMax]*1.2);
    