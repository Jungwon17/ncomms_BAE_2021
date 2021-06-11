clc; clearvars; close all; 

%% Load cell list
target = [tag.wsrxfp tag.wsefr tag.nsrxfp tag.nsefr] & T.firingRate>=0.5;
typeList = {'IT';'PT';'nsIT';'nsPT'};
nT = size(target,2);

resolution = 10;
ernum = 2;

%% Load cell data
preext = '.mat';
curext = '.t';
tFile = cellfun(@(x) regexprep(x,preext,curext), T.cellList, 'UniformOutput', false);

[pvals,selec_sxp] = deal(cell(nT,2));
selec = cell(nT,1);

for iT = 1:2
    nC = sum(target(:,iT));
    matList = T.cellList(target(:,iT));
    [tData, tList] = tLoad(tFile(target(:,iT)));
    
    fr = T.firingRate(target(:,iT));
    hyperLoc = T.hyperLocation(target(:,iT));
    
    [pvals{iT,1},pvals{iT,2},selec{iT}] = deal(NaN(nC,3));
    
    for iC = 1:nC 
        disp([typeList{iT},' cell ', num2str(iC), ' / ', num2str(nC)]);
        load([fileparts(tList{iC}), '\Events.mat'],...
            'eventTime','rewardLickTime','reward','nTrial','choiceLickTime','sample','choice');
        
        lickTime = [rewardLickTime(:,1),choiceLickTime];
        out = sum(isnan(lickTime),2)>0;
        
        trialInd = [sample==1&choice==1, sample==1&choice==2,...
            sample==2&choice==2, sample==2&choice==1];
        
        spklick = NaN(nTrial,2);
        for iPhase = 1:2
            spklick(:,iPhase) = cellfun(@length,spikeWin(tData{iC},lickTime(:,iPhase),[0 1000]));
        end
        spkrw = cellfun(@length,spikeWin(tData{iC},lickTime(:,2),[0 1500]));
        
        if sum(spklick(:))>0
            in2way = (trialInd(:,1) | trialInd(:,3))&~out;
            tbl = simple_mixed_anova(spklick(in2way,:),sample(in2way),{'phase'},{'sample'});
            pvals{iT,1}(iC,:) = tbl.pValue([2 4 5]); %sample, phase, sample*phase
        end
        
        if sum(sum(trialInd)>=ernum)==4
            outr = isnan(choiceLickTime);
            pvals{iT,2}(iC,:) = anovan(spkrw(~outr),[sample(~outr),reward(~outr)],'model','interaction','display','off');
        end
        
        spkRave = nanmean(spklick(trialInd(:,1)&~out,:),2);
        spkLave = nanmean(spklick(trialInd(:,3)&~out,:),2);
        
        spkR = spklick(trialInd(:,1)&~out,:);
        spkL = spklick(trialInd(:,3)&~out,:);
        
        if hyperLoc(iC)=='R'
            selec{iT}(iC,1) = (nanmean(spkRave)-nanmean(spkLave))/sqrt(nanstd(spkRave)^2+nanstd(spkLave)^2);
            for i = 1:2
                if isempty(selec_sxp{iT,i})
                    selec_sxp{iT,i} = NaN(nC,1);
                end
                selec_sxp{iT,i}(iC) = (nanmean(spkR(:,i))-nanmean(spkL(:,i)))/sqrt(nanstd(spkR(:,i))^2+nanstd(spkL(:,i))^2);
            end
        else
            selec{iT}(iC,1) = (nanmean(spkLave)-nanmean(spkRave))/sqrt(nanstd(spkRave)^2+nanstd(spkLave)^2);
            for i = 1:2
                if isempty(selec_sxp{iT,i})
                    selec_sxp{iT,i} = NaN(nC,1);
                end
                selec_sxp{iT,i}(iC) = (nanmean(spkL(:,i))-nanmean(spkR(:,i)))/sqrt(nanstd(spkL(:,i))^2+nanstd(spkR(:,i))^2);
            end
        end
        
        spkS = spklick((trialInd(:,1)|trialInd(:,3))&~out,1);
        spkC = spklick((trialInd(:,1)|trialInd(:,3))&~out,2);
        selec{iT}(iC,2) = (nanmean(spkS)-nanmean(spkC))/sqrt(nanstd(spkS)^2+nanstd(spkC)^2);
        
        if sum(sum(trialInd)>=ernum)==4
        spkRw = spkrw(trialInd(:,1)|trialInd(:,3));
        spkNRw = spkrw(trialInd(:,2)|trialInd(:,4));
        selec{iT}(iC,3) = (nanmean(spkRw)-nanmean(spkNRw))/sqrt(nanstd(spkRw)^2+nanstd(spkNRw)^2);
        end
    end
end


fon = cellfun(@(x) nansum(x<0.05)./sum(~isnan(x)),pvals(1:2,1),'UniformOutput',false);

nsig = cellfun(@(x) nansum(x<0.05),pvals(1:2,1),'UniformOutput',false);
nnsig = cellfun(@(x) nansum(x>0.05),pvals(1:2,1),'UniformOutput',false);
for i = 1:3
[~,pfisher(i)] = fishertest([nsig{1}(i),nnsig{1}(i);nsig{2}(i),nnsig{2}(i)]);
end

%% figure
clr = {[0.745 0.102 0.125];[0.055 0.451 0.725];[0.745 0.102 0.125];[0.055 0.451 0.725]};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 4]);
axes('Position',axpt(10,10,2:10,1:8));
hold on;
for iT = 1:2
plot(1:3,fon{iT}(1:3),'Color',clr{iT},'LineWidth',1,'Marker','o','MarkerFaceColor',[1 1 1],'MarkerSize',3)
for i = 1:3
    if pfisher(i)<0.05
        text(i,0.43,'*','FontSize',7);        
    end
end
end
xlim([0.5 3.5]);
ylim([0 0.5]);
set(gca,'XTick',1:3,'XTickLabel',{'Target';'Phase';'Interaction'},'FontSize',8,...
    'LineWidth',0.35,'YTick',0:0.1:0.5,'TickDir','out','Box','off','YTickLabel',0:10:50,...
    'XTickLabelRotation',45);
ylabel('FON (%)','FontSize',8);

fon2 = cellfun(@(x) nansum(x<0.05)./sum(~isnan(x)),pvals(1:2,2),'UniformOutput',false);

nsig2 = cellfun(@(x) nansum(x<0.05),pvals(1:2,2),'UniformOutput',false);
nnsig2 = cellfun(@(x) nansum(x>0.05),pvals(1:2,2),'UniformOutput',false);
for i = 1:3
[~,pfisher2(i)] = fishertest([nsig2{1}(i),nnsig2{1}(i);nsig2{2}(i),nnsig2{2}(i)]);
end

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 4]);
axes('Position',axpt(10,10,2:10,1:8));
hold on;
for iT = 1:2
plot(1:3,fon2{iT}(1:3),'Color',clr{iT},'LineWidth',1,'Marker','o','MarkerFaceColor',[1 1 1],'MarkerSize',3)
for i = 1:3
    if pfisher2(i)<0.05
        text(i,0.43,'*','FontSize',7);        
    end
end
end
xlim([0.5 3.5]);
ylim([0 0.5]);
set(gca,'XTick',1:3,'XTickLabel',{'Target';'Reward';'interaction'},'FontSize',8,...
    'LineWidth',0.35,'YTick',0:0.1:0.5,'TickDir','out','Box','off','YTickLabel',0:10:50,...
    'XTickLabelRotation',45);

ylabel('FON (%)','FontSize',8);
%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 4]);
axes('Position',axpt(10,10,2:10,1:8));
hold on;
for iT = 1:2
    for i = 1:2
        scatter(iT+rand(1,nnsig{iT}(i))*0.6+(i-1)*2.5,selec{iT}(pvals{iT}(:,i)>0.05,i),2,'o',...
            'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',clr{iT},'LineWidth',0.35);
        scatter(iT+rand(1,nsig{iT}(i))*0.6+(i-1)*2.5,selec{iT}(pvals{iT}(:,i)<0.05,i),6,'o',...
            'MarkerFaceColor',clr{iT},'MarkerEdgeColor','none','LineWidth',0.35);
        errorbar(iT+0.3+(i-1)*2.5,nanmean(selec{iT}(:,i)),nanstd(selec{iT}(:,i))/...
            sqrt(sum(~isnan(selec{iT}(:,i)))),'Color','k','LineWidth',0.5);
        [~,pttest(iT,i),~,stat] = ttest(selec{iT}(:,i));
        df(iT,i) = stat.df;
        tstat(iT,i) = stat.tstat;
        psignrank(iT,i) = signrank(selec{iT}(:,i));
    end
end
plot([0 10],[0 0],'k:','LineWidth',0.35);
xlim([0.7 5.4]);
ylim([-1.2 1.2]);
set(gca,'XTick',1.8:2.5:4.3,'XTickLabel',{'Sample';'Phase'},'FontSize',8,...
    'LineWidth',0.35,'YTick',-1.2:0.6:1.2,'TickDir','out','Box','off',...
    'XTickLabelRotation',45);
ylabel('Selectivity index','FontSize',8);

%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 1.5 4]);
axes('Position',axpt(10,10,2:10,1:8));
hold on;
for iT = 1:2
        scatter(iT+rand(1,sum(pvals{iT,2}(:,2)>0.05))*0.6,selec{iT}(pvals{iT,2}(:,2)>0.05,3),2,'o',...
            'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',clr{iT},'LineWidth',0.35);
        scatter(iT+rand(1,sum(pvals{iT,2}(:,2)<0.05))*0.6,selec{iT}(pvals{iT,2}(:,2)<0.05,3),6,'o',...
            'MarkerFaceColor',clr{iT},'MarkerEdgeColor','none','LineWidth',0.35);
        errorbar(iT+0.3,nanmean(selec{iT}(:,3)),nanstd(selec{iT}(:,3))/...
            sqrt(sum(~isnan(selec{iT}(:,3)))),'Color','k','LineWidth',0.75);
        [~,pttest(iT,3),~,stat] = ttest(selec{iT}(:,3));
        psignrank(iT,3) = signrank(selec{iT}(:,3));
        df(iT,3) = stat.df;
        tstat(iT,3) = stat.tstat;
        if pttest(iT,3)<0.05
           text(iT+0.3,1,'*','FontSize',8);
        end            
end
plot([0 10],[0 0],'k:','LineWidth',0.35);
xlim([0.7 2.9]);
ylim([-1.2 1.2]);
set(gca,'XTick',1.8,'XTickLabel',{'Reward'},'FontSize',8,...
    'LineWidth',0.35,'YTick',-1.2:0.6:1.2,'TickDir','out','Box','off',...
    'XTickLabelRotation',45);
ylabel('Selectivity index','FontSize',8);
