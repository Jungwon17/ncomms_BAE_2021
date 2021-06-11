function crossCorr()
% tFile = {'D:\Cheetah_data\classical_conditioning\CC-SOM-ChR3\2016-04-01_AD_1.40DV'};
close all;
[tData, tList] = tLoad();
if isempty(tList); return; end;
nT = length(tList);

load('Events.mat', 'taskTime');
tData = cellfun(@(x) x(x<taskTime(2)), tData, 'UniformOutput', false);

binSize = 1;
winWidth = 40;
gapS = [0.015 0.04];

nSubFig = 8;
nFig = ceil(nT/nSubFig);
hFig = zeros(nFig);
for iX = 1:nFig
    for iY = 1:nFig
        hFig(iX, iY) = figure;
    end
end

for iX = 1:nT
    for iY = 1:nT
        [cData, cTime] = CrossCorr(tData{iX}, tData{iY}, binSize, floor(2*winWidth/binSize)+1);
        if iX == iY
            cData(cTime==0) = 0;
        end
        
        figure(hFig(ceil(iX/nSubFig), ceil(iY/nSubFig)));
        axes('Position', axpt(nSubFig, nSubFig, mod(iX-1,nSubFig)+1, mod(iY-1,nSubFig)+1, [], gapS));
        hold on;
        hB = bar(cTime, cData, 'histc');
        set(hB, 'FaceColor', 'k', 'LineStyle', 'none');
        
        if iX ~= iY
            cJit = CrossCorrJitter(tData{iX}, tData{iY}, binSize, floor(2*winWidth/binSize)+1, 5, 1000);
            cJit = sort(cJit');
        
            lowJitter = cJit(11,:);
            highJitter = cJit(1000-10, :);
            
            plot(cTime, lowJitter, 'b-', 'LineWidth', 0.5);
            plot(cTime, highJitter, 'r-', 'LineWidth', 0.5);
        end
        
        yM = max(cData)*1.25 + 10^-10;
        title([num2str(iX),'¡æ', num2str(iY)], 'FontSize', 3);
        set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 3, 'LineWIdth', 0.2, ...
            'XLim', [-winWidth winWidth], 'XTick', [-winWidth, -5, 0, 5,winWidth], ...
            'YLim', [0 yM], 'YTick', yM);
    end
end

[cellcd,~,~] = fileparts(tList{1});
cell_filename = regexp(cellcd,'\','split');

for jX = 1:nFig
    for jY = 1:nFig
        cellfile = strcat(cell_filename(end-1),'_',cell_filename(end),'_CrossCorrelation_',num2str(jX),'_',num2str(jY),'.tif');
        print(hFig(jX,jY), '-dtiff', '-r600', cellfile{1});
    end
end