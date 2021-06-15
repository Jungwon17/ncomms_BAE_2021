clc; clearvars; close all; 
%%
% Code written by Bae JW, Jeong H 
% 2021 Nat. Commun. [Parallel processing of working memory and temporal
% information by distinct types of cortical projection neurons]
% last edited by Bae JW 2021-06-11
%%
preext = '.mat';
curext = '.t';
tFile = cellfun(@(x) regexprep(x,preext,curext), T.cellList, 'UniformOutput', false);


%% Load cell list
target = [tag.wsrxfp tag.wsefr tag.nsrxfp tag.nsefr]&T.firingRate>=0.5;
typeList = {'IT';'PT';'nsIT';'nsPT'};

[pcorr,rcorr] = deal(cell(size(target,2),1));
for iT = 1:size(target,2)
    nC = sum(target(:,iT));
    [pcorr{iT},rcorr{iT}] = deal(NaN(nC,3));
    [tData, tList] = tLoad(tFile(target(:,iT)));
    
    hyperLoc = T.hyperLocation(target(:,iT));
    for iC = 1:nC % will be replaced with for-loop
        disp([typeList{iT},' cell ', num2str(iC), ' / ', num2str(nC)]);
        load([fileparts(tList{iC}), '\Events.mat'], ...
            'taskTime','eventTime','lickOnsetTime');
        
        licktime = cat(1,lickOnsetTime{:});
        lickid = [ones(size(lickOnsetTime{1}));ones(size(lickOnsetTime{2}))*2];
        [licktime,sortIdx] = sort(licktime);
        lickid = lickid(sortIdx);
        
        lickid(licktime<=taskTime(1)) = [];
        lickid(licktime>=taskTime(2)) = [];
        licktime(licktime<=taskTime(1)) = [];
        licktime(licktime>=taskTime(2)) = [];        
        
        lickbin = [];
        lickbin(1,:) = histcounts(licktime,taskTime(1):1000:taskTime(2));
        if hyperLoc(iC)=='R'
            lickbin(2,:) = histcounts(licktime(lickid==1),taskTime(1):1000:taskTime(2));
            lickbin(3,:) = histcounts(licktime(lickid==2),taskTime(1):1000:taskTime(2));
        else
            lickbin(3,:) = histcounts(licktime(lickid==1),taskTime(1):1000:taskTime(2));
            lickbin(2,:) = histcounts(licktime(lickid==2),taskTime(1):1000:taskTime(2));
        end
        
        spkbin = histcounts(tData{iC},taskTime(1):1000:taskTime(2));
        
        for i = 1:3 %total, ipsi, contra
            [rtmp,ptmp] = corrcoef(lickbin(i,:),spkbin);
            pcorr{iT}(iC,i) = ptmp(1,2);
            rcorr{iT}(iC,i) = rtmp(1,2);
        end
    end
end


iL = 1;
clr = {[0.745 0.102 0.125];[0.055 0.451 0.725];[0.745 0.102 0.125];[0.055 0.451 0.725]};

fon = cellfun(@(x) nanmean(x<0.05),pcorr,'UniformOutput',false);

%% figure
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 1.5 4]); 
axes('Position',axpt(10,10,2:10,1:8));
hold on;
for iT = 1:2
plot(1,fon{iT}(iL),'Color',clr{iT},'LineWidth',1,'Marker','o','MarkerFaceColor',[1 1 1],'MarkerSize',3)
end

xlim([0.5 1.5])
ylim([0 1]);
set(gca,'XTick',1:3,'XTickLabel',{'sample';'phase';'interaction'},'FontSize',8,...
    'LineWidth',0.35,'YTick',0:0.25:1,'TickDir','out','Box','off','YTickLabel',0:25:100,...
    'XTickLabelRotation',45);
ylabel('FON (%)','FontSize',8);

nnsig = cellfun(@(x) sum(x>0.05),pcorr,'UniformOutput',false);
nsig = cellfun(@(x) sum(x<0.05),pcorr,'UniformOutput',false);

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 1.5 4]);
axes('Position',axpt(10,10,2:10,1:8));
hold on;
for iT = 1:2
        scatter(iT+rand(1,nnsig{iT}(iL))*0.6,rcorr{iT}(pcorr{iT}(:,iL)>0.05,iL),2,'o',...
            'MarkerFaceColor','none','MarkerEdgeColor',clr{iT},'LineWidth',0.35);
        scatter(iT+rand(1,nsig{iT}(iL))*0.6,rcorr{iT}(pcorr{iT}(:,iL)<0.05,iL),6,'o',...
            'MarkerFaceColor',clr{iT},'MarkerEdgeColor','none','LineWidth',0.35);
        errorbar(iT+0.3,nanmean(rcorr{iT}(:,iL)),nanstd(rcorr{iT}(:,iL))/...
            sqrt(sum(~isnan(rcorr{iT}(:,iL)))),'Color','k','LineWidth',0.75);
        [~,pttest(iT),~,stat] = ttest(rcorr{iT}(:,iL));
        df(iT) = stat.df;
        tstat(iT) = stat.tstat;
%         psignrank(iT) = signrank(rcorr{iT}(pcorr{iT}(:,iL)<0.05,iL));
    if pttest(iT)<0.05
        text(iT+0.3,0.55,'*','FontSize',8);
    end
end
plot([0 10],[0 0],'k:','LineWidth',0.35);
xlim([0.7 2.9]);
ylim([-0.6 0.6]);
set(gca,'XTick',[1.3 2.3],'XTickLabel',{'IT'; 'PT'},'FontSize',8,...
    'LineWidth',0.35,'YTick',-1.2:0.6:1.2,'TickDir','out','Box','off',...
    'XTickLabelRotation',45);
ylabel({'Correlation coefficient'; '(lick rate vs. firing rate)'},'FontSize',8);


ct = cbrewer('qual','Dark2',8);
ct = ct(8,:);

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6.5 3.15]);
for iT = 1:2
    axes('Position',axpt(2,1,iT,1,axpt(10,10,2:10,1:9),[0.05 0.1]))
    hold on;
    scatter(rcorr{iT}(:,2),rcorr{iT}(:,3),6,ct,'filled');
    [beta,~,stat] = glmfit(rcorr{iT}(:,2),rcorr{iT}(:,3)); 
    plot([-0.6 0.6],[-0.6 0.6]*beta(2)+beta(1),'Color','r','LineWidth',1);
    text(-0.55,0.55,['\beta = ',num2str(round(beta(2)*100)/100),'; p = ',...
        num2str(round(stat.p(2)*10000)/10000)],'FontSize',6);
    plot([0 0 NaN -1 1],[-1 1 NaN 0 0],'Color',[0.8 0.8 0.8],'LineStyle',':','LineWidth',0.35);
    xlim([-0.6 0.6])
    ylim([-0.6 0.6]);
    set(gca,'XTick',-0.6:0.6:0.6,'YTick',-0.6:0.6:0.6,'FontSize',8,'Box','off','TickDir','out','LineWidth',0.35);
     
    title(typeList{iT},'FontSize',8,'FontWeight','Bold','Color',clr{iT});
    if iT==1
        ylabel({'Ipsilateral';'Correlation coefficient'},'FontSize',8);
        xlabel({'Correlation coefficient';'Contralateral'},'FontSize',8)
    else
       set(gca,'YTickLabel',[]);
    end
end