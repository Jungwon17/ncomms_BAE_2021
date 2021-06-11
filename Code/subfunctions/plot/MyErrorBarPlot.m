function h = MyScatterBarPlot(y,x,barWidth, xColor, sigGroup, yLim)
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
group = unique(x);
nGroup = length(group);

hold on;
yMax = zeros(1, nGroup);
yMin = zeros(1, nGroup);
for iGroup = 1:nGroup
    yPoint = y(x==group(iGroup) & ~isnan(y));
    nPoint = sum(~isnan(yPoint));
    xPoint = iGroup + 0.75*barWidth*(rand(nPoint,1)-0.5);
    
    yMean = nansum(yPoint)/nPoint;
    ySem = nanstd(yPoint)/sqrt(nPoint);
    yMax(iGroup) = yMean + ySem;
    yMin(iGroup) = yMean - ySem;
    
    h(iGroup).bar = bar(iGroup,yMean,'FaceColor',xColor{iGroup},'LineStyle','none','BarWidth',barWidth);
    h(iGroup).errorbar = errorbar(iGroup,yMean,ySem,'LineWidth',2,'Color',xColor{iGroup});
    
    errorbarT(h(iGroup).errorbar,0.25,0.5);
end

if nargin < 6
    if min(yMin)>0
        yLim = [0 ceil(max(yMax)*1.2)];
    else
        yLim = [floor(min(yMin)*1.2) ceil(max(yMax)*1.2)];
    end
end

if nargin>=5 && iscell(sigGroup) && ~isempty(sigGroup)
    nG = length(sigGroup);
    for iG = 1:nG
        plot([sigGroup{iG}(1) sigGroup{iG}(1)], [yMax(sigGroup{iG}(1))+diff(yLim)*0.025 max(yMax(sigGroup{iG}))+diff(yLim)*0.05], ...
            'LineWidth', 0.35, 'Color', 'k');
        plot([sigGroup{iG}(2) sigGroup{iG}(2)], [yMax(sigGroup{iG}(2))+diff(yLim)*0.025 max(yMax(sigGroup{iG}))+diff(yLim)*0.05], ...
            'LineWidth', 0.35, 'Color', 'k');
        plot([sigGroup{iG}(1) sigGroup{iG}(2)], [max(yMax(sigGroup{iG}))+diff(yLim)*0.05 max(yMax(sigGroup{iG}))+diff(yLim)*0.05], ...
            'LineWidth', 0.35, 'Color', 'k');
        text(mean(sigGroup{iG}), max(yMax(sigGroup{iG}))+diff(yLim)*0.075, '*', ...
            'FontSize', 4, 'HorizontalAlignment', 'center');
    end
end

set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.2,...
    'XLim',[0.5 nGroup+0.5],'XTick',1:nGroup, ...
    'YLim', yLim, 'YTick', unique([0 yLim]));
