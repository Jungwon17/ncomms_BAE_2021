close all; clc; clearvars;
%%
% Temporal decoding with SVM method

% Code written by Bae JW, Jeong H 
% 2021 Nat. Commun. [Parallel processing of working memory and temporal
% information by distinct types of cortical projection neurons]
% last edited by Bae JW 2021-06-11

Directory='D:\Backup\code\WM_2021';
cd(Directory); % set your directory here
%% load cell data

cd([Directory,'\Data\IT'])
IT_Neurons = FindFiles('T*.mat','CheckSubdirs',1);
[lvcd,name] = cellfun(@fileparts,IT_Neurons,'UniformOutput',false);

cd([Directory,'\Data\PT'])
PT_Neurons = FindFiles('T*.mat','CheckSubdirs',1);
[lvcd2,name2] = cellfun(@fileparts,IT_Neurons,'UniformOutput',false);
%%
nT = 2;
nEnsemble = 26;
win = [2500 6000];
binSize = (win(2)-win(1))/10; 
nBin = diff(win)/binSize;
nTrain = 9;
nTest = 1;
nIter = 100;
typeList = {'IT';'PT'};
[result,accuracy,r,p,errorbin,P_corr] = deal(cell(nT,1));

for iT = 1:nT
    
    if iT == 1 % IT
        cellList = IT_Neurons;
        cd([Directory,'\Data\IT']);
        load('hyperLoc');
    else
        cellList = PT_Neurons;
        cd([Directory,'\Data\PT']);
        load('hyperLoc');
    end

    nC = length(cellList);
    spkTotal = cell(nC,nBin);
    for iC = 1:nC
        load(cellList{iC},'spikeTime');
        load([fileparts(cellList{iC}),'\Events.mat'],'trialIndex','trialResult');
        
        [time,spktmp] = spikeBin(spikeTime(trialIndex(:,1)|trialIndex(:,4)),win,binSize,binSize);
        spkTotal(iC,:) = mat2cell(spktmp,trialResult(1)+trialResult(4),ones(1,length(time)));
    end
    group = repmat(time,nTrain,1);
    [result{iT},accuracy{iT}] = deal(NaN(nIter,nTrain+nTest,length(time)));
    [r{iT},p{iT}] = deal(NaN(nIter,1));
    errorbin{iT} = NaN(nIter,nBin);
    P_corr{iT} = NaN(nIter,nBin);
    for iIter = 1:nIter
        ensembleInd = randsample(nC,nEnsemble);
        spkIter = cellfun(@(x) x(randsample(size(x,1),nTrain+nTest)),spkTotal(ensembleInd,:),...
            'UniformOutput',false)';
        error = NaN(nTrain+nTest,nBin);
        for jIter = 1:(nTrain+nTest)
            fprintf('%s: %d/%d subiteration of %d/%d iteration\n',...
                typeList{iT},jIter,nTrain+nTest,iIter,nIter);
            if jIter==1
                trainInd = 2:(nTrain+nTest);
            elseif jIter==nTrain+nTest
                trainInd = 1:nTrain+nTest-1;
            else
                trainInd = [1:jIter-1, jIter+1:(nTrain+nTest)];
            end
            spkTrain = cell2mat(cellfun(@(x) x(trainInd),spkIter,'UniformOutput',false));
            spkTest = cell2mat(cellfun(@(x) x(jIter),spkIter,'UniformOutput',false));
            mdl = fitcecoc(spkTrain,group(:));
            result{iT}(iIter,jIter,:) = predict(mdl,spkTest);
            accuracy{iT}(iIter,jIter,:) = predict(mdl,spkTest)==time';
            error(jIter,:) = abs(predict(mdl,spkTest)-time');
        end
        data = squeeze(result{iT}(iIter,:,:));
        g = repmat(time,nTrain+nTest,1);
        [rtmp,ptmp] = corrcoef(data(:),g(:));
        r{iT}(iIter) = rtmp(1,2);
        p{iT}(iIter) = ptmp(1,2);
        errorbin{iT}(iIter,:) = nanmean(error)/binSize;
        P_corr{iT}(iIter,:)=sum(data(:)==g(:))/size(data(:),1);
        
    end
end


%% figure
win = [2500 6000];
binSize = (win(2)-win(1))/10; 
nBin = diff(win)/binSize;
nIter = 100;

clr = {[0.745 0.102 0.125];[0.055 0.451 0.725]};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 4]);
typeList = {'IT';'PT'};
% typeList = {'FS';'PP'};
for iT = 1:2
    match{iT} = zeros(nIter,nBin,nBin);
    for iIter = 1:100
        for jIter = 1:20
            for iTime = 1:nBin
                match{iT}(iIter,find(result{iT}(iIter,jIter,iTime)==time),iTime) =...
                    match{iT}(iIter,find(result{iT}(iIter,jIter,iTime)==time),iTime)+1;
            end
        end
        maxmatch = max(squeeze(match{iT}(iIter,:,:)));
        match{iT}(iIter,:,:) = squeeze(match{iT}(iIter,:,:))./repmat(maxmatch,nBin,1); % winner takes all, max normalized
%         match{iT}(iIter,:,:) = match{iT}(iIter,:,:)/(nTrain+nTest);
    end

    axes('Position',axpt(2,1,iT,1,axpt(1,1,1,1,axpt(15,15,2:14,1.5:12.5)),[0.05 0.05]))
    imagesc(time,time,squeeze(nanmean(match{iT},1)));
    axis xy
    
    set(gca,'CLim',[0 1],'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,...
        'XTick',[2500,6000],'XTickLabel',[0.5 4],'YTick',[2500,6000],'YTickLabel',[0.5 4]);
    
    if iT~=1
        set(gca,'YTickLabel',[]);
    else
        ylabel({'Predicted';'Time from delay onset (s)'},'FontSize',5);
    end
    xlabel({'Time from delay onset (s)';'Actual'},'FontSize',5);
    title(typeList{iT},'FontSize',7,'Color',clr{iT});
end

axes('Position',axpt(15,1,15,1))
c = colorbar('Ticks',[0 0.5 1],'TickLength',0.01,'YAxisLocation','right'...
    ,'AxisLocation','out','TickDirection','out','Box','off','FontSize',5);
set(gca,'visible','off','CLim',[0 1]);
c.Position = [0.91 0.5825 0.03 0.3];
c.FontSize = 5;
c.Box ='off';
axis off
text(c.Position(1)-2,c.Position(2)+0.42,{'Readout';'value'},'FontSize',4)
%   Read out value : 10.1523/JNEUROSCI.1789-16.2016  

%% fig

clr = {[0.745 0.102 0.125];[0.055 0.451 0.725];[0.2 0.2 0.2];[0.7 0.7 0.7]};

% clr = {[0.7 0.7 0.7];[0.2 0.2 0.2]};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 4]);
axes('Position',axpt(1,1,1,1,axpt(1,1,1,1,axpt(15,15,2.5:14,1.5:12.5)),[0.05 0.05]))
hold on;
for i =1:2
    if i==1; ii=1; else ii=3; end
    si= r2{i}';
    si=si';nanmean_si=nanmean(si); si(isnan(si))=[];s21=std(si);sem21=s21/sqrt(length(si));
    bar([ii],nanmean_si,0.5,'FaceColor',clr{i},'EdgeColor',clr{i});hold on;
    alpha(0.6);
    errorbar([ii],[nanmean_si ],[sem21], '.','Color',clr{i}-[0.05 0.05 0.05],'LineWidth', 1,'CapSize',5);
    
    si= r{i}';
    si=si';nanmean_si=nanmean(si); si(isnan(si))=[];s21=std(si);sem21=s21/sqrt(length(si));
    bar([ii+0.8],nanmean_si,0.5,'FaceColor',clr{i},'EdgeColor',clr{i},...
        'Linestyle','--');hold on;
    alpha(0.6);
    errorbar([ii+0.8],[nanmean_si ],[sem21], '.','Color',clr{i}-[0.05 0.05 0.05],'LineWidth', 1,'CapSize',5);

end

plot([1 3],[0.7 0.7],'Color',[0 0 0],...
    'LineWidth',0.35);
plot([1 1.8],[0.2 0.2],'Color',[0 0 0],...
    'LineWidth',0.3);
plot([3 3.8],[0.6 0.6],'Color',[0 0 0],...
    'LineWidth',0.3);

text(1.7,0.73,'***','FontSize',8);
text(1.2,0.23,'###','FontSize',4);
text(3.2,0.63,'###','FontSize',4);

xlim([0.5 4.3])
ylim([-0.05 0.8])

set(gca,'Box','off','TickDir','out','LineWidth',0.35,'FontSize',5,...
    'XTick',[1, 1.8 3 3.8],'XTickLabel',{'IT';'Shifted-IT';'PT';'Shifted-PT'},'YTick',[0 0.4 0.8]);
ylabel('r (actual vs. predictd)','FontSize',5);
xtickangle(45)
