function h = MyErrorBarGroupPlot(y,x,barWidth,xColor, sigGroup, yLim)
%MyScatterBarPlot Plot scatter plot and Bar plot with sem error bar
%
% function MyScatterBarPlot(data,group,groupColor)
%
%
% Inputs
% y - value
% x - group variable for y
%   (first column for different groups, two column for different strata)
%
%   h = MyErrorBarGroupPlot(y, x, 0.9, fillColor, 1:4, [0.3 1.05]);
%
%   Only compatible with <2014a version
%
% Author: Dohoung Kim
% Created: June 2015
% Updated: May 2016
if nargin < 6
    sigGroup = [];
end

group = unique(x(:, 1));
strata = unique(x(:, 2));

nGroup = length(group);
nStrata = length(strata);


[yMean, ySem] = deal(zeros(nGroup, nStrata));
for iGroup = 1:nGroup
    for iStrata = 1:nStrata
        yPoint = y(x(:, 1)==group(iGroup) & x(:, 2)==strata(iStrata) & ~isnan(y));
        nPoint = sum(~isnan(yPoint));
        
        yMean(iGroup, iStrata) = nansum(yPoint)/nPoint;
        ySem(iGroup, iStrata) = nanstd(yPoint)/sqrt(nPoint);
    end  
end

hold on;
h.bar = bar(yMean, 'BarWidth', barWidth);

xpt = zeros(nGroup, nStrata);
for iStrata = 1:nStrata
    set(h.bar(iStrata), 'FaceColor', xColor{iStrata}, 'LineStyle', 'none', 'ShowBaseLine', 'off');
    
    xptTemp = get(get(h.bar(iStrata), 'Children'), 'XData');
    xpt(:, iStrata) = mean(xptTemp([1 4], :), 1);

    h.errorbar(iStrata) = errorbar(xpt(:, iStrata),yMean(:, iStrata) ,ySem(:, iStrata), ...
        'Color',xColor{iStrata});
    errorbarT(h.errorbar(iStrata), 0.1, 0.35);
end

ypt = yMean+ySem + 0.025 * diff(yLim);
for iP = 1:4
    text(xpt(iP, 1), ypt(iP, 1) + 0.025 * diff(yLim), '*', 'FontSize', 10, 'HorizontalAlignment', 'center', 'Color', [1 0 0]);
end
for iP = sigGroup
    text(xpt(iP, 2), ypt(iP, 2) + 0.025 * diff(yLim), '*', 'FontSize', 10, 'HorizontalAlignment', 'center', 'Color', [1 0 0]);
end

% yPeak = max(yMean+ySem, [], 2) + 0.05 * diff(yLim);
% ypt = yMean+ySem + 0.025 * diff(yLim);
% for iP = sigGroup
%     plot([xpt(iP, 1) xpt(iP, 1)], [ypt(iP, 1), yPeak(iP)], 'LineWidth', 0.35, 'Color', 'k');
%     plot([xpt(iP, 2) xpt(iP, 2)], [ypt(iP, 2), yPeak(iP)], 'LineWidth', 0.35, 'Color', 'k');
%     plot([xpt(iP, 1) xpt(iP, 2)], [yPeak(iP) yPeak(iP)], 'LineWidth', 0.35, 'Color', 'k');
%     text(iP, yPeak(iP)+ 0.025 * diff(yLim), '*', 'FontSize', 5, 'HorizontalAlignment', 'center');
% end

set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,...
    'XLim',[0.5 nGroup+0.5],'XTick',1:nGroup);
