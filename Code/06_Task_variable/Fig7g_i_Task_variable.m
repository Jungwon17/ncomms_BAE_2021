clc; clearvars; close all; 

%% Load cell list
target = [tag.wsrxfp tag.wsefr tag.pc] & T.firingRate>=0.5;
typeList = {'IT';'PT';'PP'};
nT = size(target,2);
resolution = 10;

preext = '.mat';
curext = '.t';
tFile = cellfun(@(x) regexprep(x,preext,curext), T.cellList, 'UniformOutput', false);

%%
[panova,fvalue,pttest,psthconv,selec] = deal(cell(2,3));
for iT = 1:2
    nC = sum(target(:,iT));
    matList = T.cellList(target(:,iT));
    [tData, tList] = tLoad(tFile(target(:,iT)));
    
    fr = T.firingRate(target(:,iT));
    hyperLoc = T.hyperLocation(target(:,iT));

    for iC = 1:nC 
        disp([typeList{iT},' cell ', num2str(iC), ' / ', num2str(nC)]);
        load([fileparts(tList{iC}), '\Events.mat'],...
            'eventTime','rewardLickTime','reward','nTrial','choiceLickTime','sample','choice');
        
        onsetTime = eventTime(:,[1,4,7]);
        out = sum(isnan(onsetTime),2)>0;
        
        for iE = 1:3
            if isempty(pttest{iT,iE})
                [pttest{iT,iE},selec{iT,iE}] = deal(NaN(nC,1));
                [panova{iT,iE},fvalue{iT,iE}] = deal(NaN(nC,3));
            end
            
           spikeTime = spikeWin(tData{iC},onsetTime(:,iE),[-2000 2000]);
           [~,spk] = spikeBin(spikeTime,[-500 500],500,500);
           [~,pttest{iT,iE}(iC)] = ttest(spk(~out,1),spk(~out,2));
           
           tbl = simple_mixed_anova(spk(~out,:),sample(~out),{'phase'},{'sample'});
           panova{iT,iE}(iC,:) = tbl.pValue([2 4 5]); %sample, phase, sample*phase
           fvalue{iT,iE}(iC,:) = tbl.F([2 4 5]);
           
           selec{iT,iE}(iC) = (nanmean(spk(~out,2))-nanmean(spk(~out,1)))/...
               sqrt(nanstd(spk(~out,2))^2+nanstd(spk(~out,1))^2);
           
           [time,spk] = spikeBin(spikeTime,[-2000 2000],10,10);
           spk = spk*(1000/10);
           if isempty(psthconv{iT,iE})
              psthconv{iT,iE} = deal(NaN(nC,3,length(time)));
           end
           psthconv{iT,iE}(iC,1,:) = conv(nanmean(spk(~out,:)),fspecial('Gaussian',[1 5*resolution],resolution),'same');
           for iSample = 1:2
               convtmp = conv(nanmean(spk(~out&sample==iSample,:)),fspecial('Gaussian',[1 5*resolution],resolution),'same');
               if hyperLoc(iC)=='R'
                   psthconv{iT,iE}(iC,iSample+1,:) = convtmp;
               else
                   psthconv{iT,iE}(iC,4-iSample,:) = convtmp;
               end
           end
        end
    end
end

nsig = cellfun(@(x) sum(x<0.05),pttest);
nnsig = cellfun(@(x) sum(x>0.05),pttest);
fon = cellfun(@(x) nanmean(x<0.05),pttest);
clr = {[0.745 0.102 0.125];[0.055 0.451 0.725];[0.6 0.6 0.6]};

%% figure
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 4]);
hold on;
for iT = 1:2
   plot(1:3,fon(iT,:),'Color',clr{iT},'Marker','o','MarkerSize',3,'MarkerFaceColor',[1 1 1],'LineWidth',1); 
   pp = chisq([cellfun(@(x) sum(x<0.05),pttest(iT,:));cellfun(@(x) sum(x>0.05),pttest(iT,:))]');
end
for i = 1:3
   [~,pchi(i)] = chisq2(sum(pttest{1,i}>0.05),sum(pttest{1,i}<0.05),...
       sum(pttest{2,i}>0.05),sum(pttest{2,i}<0.05),0.05); 
   if pchi(i)<0.05
      text(i,0.57,'*','FontSize',8); 
   end
end
xlim([0.5 3.5])
ylim([0 0.7]);
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',1:3,'XTickLabel',{'LED';'Buzzer';'ITI'},'YTick',0:0.2:0.6,...
    'YTickLabel',0:20:60,'XTickLabelRotation',45);
ylabel('FON (%)','FontSize',8);


ylabelList = {'LED';'Buzzer';'ITI'};
binRange = -1:0.1:1;
typeList = {'IT';'PT'};

nsig = cellfun(@(x) sum(x<0.05),pttest);
nnsig = cellfun(@(x) sum(x>0.05),pttest);

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 4]);
axes('Position',axpt(10,10,2:10,1:8));
hold on;
for iT = 1:2
    for i = 1:3
        scatter(iT+rand(1,nnsig(iT,i))*0.6+(i-1)*2.5,selec{iT,i}(pttest{iT,i}>0.05),2,'o',...
            'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',clr{iT},'LineWidth',0.35);
        scatter(iT+rand(1,nsig(iT,i))*0.6+(i-1)*2.5,selec{iT,i}(pttest{iT,i}<0.05),6,'o',...
            'MarkerFaceColor',clr{iT},'MarkerEdgeColor','none','LineWidth',0.35);
        errorbar(iT+0.3+(i-1)*2.5,nanmean(selec{iT,i}),nanstd(selec{iT,i})/...
            sqrt(sum(~isnan(selec{iT,i}))),'Color','k','LineWidth',0.75);
%         [~,p_ttest(iT,i)] = ttest(selec{iT,i});
        [~,p_ttest(iT,i),~,stat] = ttest(selec{iT,i});
        df(iT,i) = stat.df;
        tstat(iT,i) = stat.tstat;
        
        if p_ttest(iT,i)<0.001
            text(iT+(i-1)*2.5,1,'***','FontSize',8);
        elseif p_ttest(iT,i)<0.01
            text(iT+0.2+(i-1)*2.5,1,'**','FontSize',8);
        elseif p_ttest(iT,i)<0.05
            text(iT+0.3+(i-1)*2.5,1,'*','FontSize',8);
%         elseif p_signrank(iT,i)<0.1
%             text(iT+0.3+(i-1)*2.5,0.9,'~','FontSize',8);
        end
    end
end
plot([0 10],[0 0],'k:','LineWidth',0.35);
xlim([0.8 7.8]);
ylim([-0.65 1.1]);
set(gca,'XTick',1.8:2.5:7.3,'XTickLabel',{'LED';'Buzzer';'ITI'},'FontSize',8,...
    'LineWidth',0.35,'YTick',-0.5:0.5:1,'TickDir','out','Box','off',...
    'XTickLabelRotation',45);
ylabel('Selectivity index','FontSize',8);


pair = [1 2; 1 3; 2 3];
pairList = {'L vs. B';'L vs. P';'B vs. P'};

stl = {'-';'--';':'};
ct = cbrewer('qual','Dark2',8);
ct = ct(8:-1:6,:);

clr2 = {'r';'k';'k'};

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6.5 3.15]);
for iT = 1:2
    axes('Position',axpt(2,1,iT,1,axpt(10,10,2:10,1:9),[0.05 0.1]))
    hold on;
    for ipair = 1:3
        if ipair==1
            scatter(selec{iT,pair(ipair,1)},selec{iT,pair(ipair,2)},6,ct(1,:),'filled');
%             p = sum(cell2mat(psignrank(iT,pair(ipair,:)))<0.05,2);
%         scatter(selec{iT,pair(ipair,1)}(p==0),selec{iT,pair(ipair,2)}(p==0),2,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.6 0.6 0.6]);
%         scatter(selec{iT,pair(ipair,1)}(p==1),selec{iT,pair(ipair,2)}(p==1),6,[0.6 0.6 0.6],'filled');
%         scatter(selec{iT,pair(ipair,1)}(p==2),selec{iT,pair(ipair,2)}(p==2),6,'k','filled');
        end
        [beta,~,stat] = glmfit(selec{iT,pair(ipair,1)},selec{iT,pair(ipair,2)});
%         if ipair==1
        plot([-1 1],[-1 1]*beta(2)+beta(1),'Color',clr2{ipair},'LineStyle',stl{ipair});
%         else
%             plot([-1 1],[-1 1]*beta(2)+beta(1),'Color',ct(ipair,:));
%         end
        plot([-0.75 -0.65],repmat(0.95-(ipair-1)*0.15,1,2),'Color','k','LineStyle',stl{ipair});
        text(-0.6,0.95-(ipair-1)*0.15,[pairList{ipair},': \beta = ',num2str(round(beta(2)*100)/100)],'FontSize',5);
    end
    plot([0 0 NaN -1 1],[-1 1 NaN 0 0],'Color',[0.8 0.8 0.8],'LineStyle',':','LineWidth',0.35);
    xlim([-0.8 1])
    ylim([-0.8 1]);
    set(gca,'XTick',-0.8:0.8:1,'YTick',-0.8:0.8:1,'FontSize',8,'Box','off','TickDir','out','LineWidth',0.35);
     xlabel({'Selectivity index'; 'LED'},'FontSize',8)
    title(typeList{iT},'FontSize',8,'FontWeight','Bold','Color',clr{iT});
    if iT==1
        ylabel({'Buzzer';'Selectivity index'},'FontSize',8);
        
    else
       set(gca,'YTickLabel',[]);
    end
end

eventList = {'LED';'Buzzer';'Port-out'};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 9.5 3.15]);
for ipair = 1:3
    axes('Position',axpt(3,1,ipair,1,axpt(1,10,1,1:9),[0.05 0.1]))
    hold on;
    for iT = 1:2
            scatter(selec{iT,pair(ipair,1)},selec{iT,pair(ipair,2)},6,clr{iT},'filled');
        [beta,~,stat] = glmfit(selec{iT,pair(ipair,1)},selec{iT,pair(ipair,2)});
        plot([-1 1],[-1 1]*beta(2)+beta(1),'Color',clr{iT},'LineWidth',1);
%         plot([-0.95 -0.8],repmat(0.95-(ipair-1)*0.15,1,2),'k','LineStyle',stl{ipair});
        text(-0.75,0.95-(iT-1)*0.15,['\beta = ',num2str(round(beta(2)*100)/100),...
            '; p = ',num2str(round(stat.p(2)*10000)/10000)],'FontSize',5,'Color',clr{iT});
    end
    if ipair==2
       title('Selectivity index','FontSize',8); 
    end
    plot([0 0 NaN -1 1],[-1 1 NaN 0 0],'Color','k','LineStyle',':','LineWidth',0.35);
    xlim([-0.8 1])
    ylim([-0.8 1]);
    set(gca,'XTick',-0.8:0.8:1,'YTick',-0.8:0.8:1,'FontSize',8,'Box','off','TickDir','out','LineWidth',0.35);
     xlabel(eventList{pair(ipair,1)},'FontSize',8)
     ylabel(eventList{pair(ipair,2)},'FontSize',8);
    if iT>1
       set(gca,'YTickLabel',[]);
    end
end