clearvars; clc; close all;

%%
% This code is modified from the code from https://github.com/machenslab/dPCA
% you need sub-functions from https://github.com/machenslab/dPCA
% last edited by JWBAE 2021-06-11

load('D:\Backup\code\WM_2021\Data\Untagged\Untagged_WS.mat');

%% dpca

combinedParams = {{1, [1 2]}, {2}};
margNames = {'Sample';'Time'};
margColours = [[0 0 0];[1 0 0]];
timeEvents = [0 2000 6000];

[W,V,whichMarg] = dpca(spkAveTotal(:,1:2,:), 20, ...
    'combinedParams', combinedParams);
explVar = dpca_explainedVariance(spkAveTotal(:,1:2,:), W, V, ...
    'combinedParams', combinedParams);
dpca_plot(spkAveTotal(:,1:2,:), W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);


%% Figure 

margNames = {'Sample-dep.';'Sample-indep.'};
clr = {[0 0 0];[0.5 0.5 0.5]};
c = [0.2422 0.1504 0.6603; 0.9769 0.9839 0.0805];
c1 = [0.8422 0.7504 1; 1 1 0.6805];
ylimit = [-0.4 0.5;-0.15 0.1;-0.4 0.4;-0.4 0.4;-0.6 0.5];
for iT = 1:nT
    fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 14 5]);
    z = cell(2,1);
    z{1} = bsxfun(@minus,squeeze(spkAveTotal{iT}(:,1,:)),mean(squeeze(spkAveTotal{iT}(:,:)),2))'*W{iT};
    z{2} = bsxfun(@minus,squeeze(spkAveTotal{iT}(:,2,:)),mean(squeeze(spkAveTotal{iT}(:,:)),2))'*W{iT};
    for iM = 1:2
        m = find(whichMarg{iT}==iM,3,'first');
        for iPC = 1:3
            axes('Position',axpt(3,2,iPC,iM,axpt(4,12,2:4,1:11,[],[0.1 0.1]),[0.05 0.1]));
            hold on;
            plot([0 0 NaN 2000 2000 NaN 6000 6000],[ylimit(iT,:), NaN ylimit(iT,:),...
                NaN ylimit(iT,:)],'Color',[0.6 0.6 0.6],'LineWidth',0.35);
            for i = 1:2
               plot(time,z{i}(:,m(iPC)),'Color',clr{i},'LineWidth',1); 
            end
            h = title(['dPC ',num2str(m(iPC)),' (',...
                num2str(round(explVar{iT}.componentVar(m(iPC))*10)/10),'%)'],...
                'FontSize',5);
            h.BackgroundColor = c1(iM,:);
            xlim([-1000 9500]);
            ylim(ylimit(iT,:));
            set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,...
                'XTick',0:2000:8000,'XTickLabel',0:2:8,'YTick',[ylimit(iT,1) 0 ylimit(iT,2)],...
                'YTickLabel',[ylimit(iT,1) 0 ylimit(iT,2)]*100);
            if iM==2
                if iPC==2
                    xlabel('Time from trial onset (s)','FontSize',5);
                end
            else
                set(gca,'XTickLabel',[]);
            end
            
            if iPC==1
                ylabel({margNames{iM};'Activity (a.u.)'},'FontSize',5);
            else
                set(gca,'YTickLabel',[]);
            end
        end
    end
    axes('Position',axpt(1,2,1,1,axpt(4,12,1,1:11,[],[0.1 0.1]),[0.05 0.1]));
    plot(1:10,cumsum(explVar{iT}.componentVar(1:10)),'ks-','MarkerSize',2);
    xlim([0 11]);
    set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,...
        'XTick',[1 10],'YTick',40:20:100,'XTickLabel',[]);
    ylabel('Explained variance (%)','FontSize',5);
    
  
    axes('Position',axpt(1,2,1,2,axpt(4,12,1,1:11,[],[0.1 0.05])));
    hbar = bar(explVar{iT}.margVar(:,1:10)','stacked');
    xlim([0 11]);
    set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,...
        'XTick',[1 10],'YTick',0:20:60);
    rectangle('Position',[4 56 1.2 2],'FaceColor',c(1,:));
    rectangle('Position',[4 51 1.2 2],'FaceColor',c(2,:));
    text(5.6,57,margNames{1},'FontSize',5);
    text(5.6,52,margNames{2},'FontSize',5);
    alpha(0.4);
    ylabel('Explained variance (%)','FontSize',5);
    xlabel('dPC','FontSize',5);
    
    axes('Position',axpt(5,4,3:5,2:3,axpt(1,2,1,2,axpt(4,12,1,1:11,[],[0.1 0.1]),[0.05 0.1])));
    p = pie(sum(explVar{iT}.margVar,2));
    t = p(2);
    t.FontSize = 5;
    t = p(4);
    t.FontSize = 5;
    alpha(0.4);
    

end


%%

z = cell(4,1);
for iTI = 1:4
    z{iTI} = NaN(nBin,20);
    z{iTI} = bsxfun(@minus,squeeze(spkAveTotal(:,iTI,:)),...
        nanmean(spkAveTotal(:,t,:),2))'*W;
end

% ylimit = {[-10 10],[-10 10];[-40 30],[-20 25]};
% ylimit = {[-5 5],[-5 5];[-10 10],[-10 10]};
ylimit = {[-10 10],[-10 10];[-30 30],[-30 30]};

c = [0 0 0; 0.5 0.5 0.5; 0 0 0; 0.5 0.5 0.5];
c1 = [0.8422 0.7504 1; 1 1 0.6805];
ylabelList = {'Sample-dep.';'Sample-indep'};

lineStl = {'-';':'};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 6]);
for iM = 1:2
    pcnumbers = find(whichMarg==iM,2,'first');
    for iPC = 1:2
        axes('Position',axpt(2,2,iPC,iM,axpt(15,15,2:15,1:14),[0.1 0.15]));
        plot([0 0 NaN 2000 2000 NaN 6000 6000],[ylimit{iM,iPC} NaN ylimit{iM,iPC} NaN ylimit{iM,iPC}],...
            'Color',[0.8 0.8 0.8],'LineWidth',0.35);
        hold on;
        for iTI = 1:4
            plot(time,z{iTI}(:,pcnumbers(iPC)),'Color',c(iTI,:),...
                'LineStyle',lineStl{ceil(iTI/2)},'LineWidth',1);
        end
        ylim(ylimit{iM,iPC});
        xlim([-1000 9500]);
        set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,...
            'XTick',0:2000:8000,'XTickLabel',0:2:8,'YTick',[ylimit{iM,iPC}(1) 0 ylimit{iM,iPC}(2)]);
        title(['dPC ',num2str(pcnumbers(iPC)),' (',...
            num2str(round(explVar.componentVar(pcnumbers(iPC))*10)/10),'%)'],...
            'FontSize',5,'BackgroundColor',c1(iM,:));    
        if iPC==1
            ylabel({ylabelList{iM};'Activity (a.u.)'},'FontSize',5);
        end
        if iM==1
            set(gca,'XTickLabel',[]);
        else
            xlabel('Time from trial onset (s)','FontSize',5);
        end
    end
end
