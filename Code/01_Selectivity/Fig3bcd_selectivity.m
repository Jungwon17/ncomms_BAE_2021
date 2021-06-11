clearvars; clc; close all;
%%
% edited from, HJ code 2020.
% last edited by JWBAE 2021-06-11

Directory='D:\Backup\code\WM_2021';
cd(Directory); % set your directory here

%% load cell data

cd([Directory,'\Data\IT'])
IT_Neurons = FindFiles('T*.mat','CheckSubdirs',1);
[lvcd,name] = cellfun(@fileparts,IT_Neurons,'UniformOutput',false);

cd([Directory,'\Data\PT'])
PT_Neurons = FindFiles('T*.mat','CheckSubdirs',1);
[lvcd2,name2] = cellfun(@fileparts,IT_Neurons,'UniformOutput',false);

%% selectivity
nT = 2;
[d,d4s,d4s_err,p] = deal(cell(nT,1));
[spkTotalDelay,spkTotalDelayz] = deal(cell(nT,2));
resolution = 10;

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
    [d4s{iT},d4s_err{iT}] = deal(NaN(nC,1));
    
    for iC = 1:nC
        iC
        load(cellList{iC},'spikeTime')
        load([fileparts(cellList{iC}),'\Events.mat'],'trialIndex','trialResult');
        
        if trialResult(1)<20 | trialResult(4)<20
            continue;
        end
        
        trialInd = trialIndex(:,[1 4]);
        trialInd_err = trialIndex(:,[2 5]);
        if hyperLoc(iC)=='L'
            trialInd = flip(trialInd,2);
            trialInd_err = flip(trialInd_err,2);
        end
        
        [time,spk] = spikeBin(spikeTime,[-1500 12000],1000,100);
        [~,spk4s] = spikeBin(spikeTime,[2000 6000],1000,1000);
        spkRef = spk4s(trialInd(:,1)|trialInd(:,2),:);
        refmean = nanmean(spkRef(:));
        refstd = nanstd(spkRef(:));
        
        spkIpsi = spk(trialInd(:,1),:);
        spkContra = spk(trialInd(:,2),:);
        if isempty(d{iT})
            [d{iT},p{iT}] = deal(NaN(nC,length(time)));
        end
        d{iT}(iC,:) = [nanmean(spkIpsi)-nanmean(spkContra)]./sqrt(nanstd(spkIpsi).^2+nanstd(spkContra).^2);
        [~,p{iT}(iC,:)] = ttest2(spkIpsi,spkContra);
        
        spkIpsi_4s = sum(spk4s(trialInd(:,1),:),2);
        spkContra_4s = sum(spk4s(trialInd(:,2),:),2);
        d4s{iT}(iC) = (nanmean(spkIpsi_4s)-nanmean(spkContra_4s))./sqrt(nanstd(spkIpsi_4s).^2+nanstd(spkContra_4s).^2);
        
        [psthtime,psth] = spikeBin(spikeTime,[0 8000],10,10);
        psth = psth*(1000/10);
        for i = 1:2
            if isempty(spkTotalDelay{iT,i})
                [spkTotalDelay{iT,i},spkTotalDelayz{iT,i}] = deal(NaN(nC,sum(psthtime>=2000 & psthtime<=6000)));
            end
            spkave = nanmean(psth(trialInd(:,i),:));
            spkconv = conv(spkave,fspecial('Gaussian',[1 5*resolution],resolution),'same');
            spkTotalDelay{iT,i}(iC,:) = spkconv(psthtime>=2000 & psthtime<=6000);
        end
        
        if sum(sum(trialInd_err)>1)==2
            spkIpsi_4s = sum(spk4s(trialInd_err(:,1),:),2);
            spkContra_4s = sum(spk4s(trialInd_err(:,2),:),2);
            d4s_err{iT}(iC) = (nanmean(spkIpsi_4s)-nanmean(spkContra_4s))./sqrt(nanstd(spkIpsi_4s).^2+nanstd(spkContra_4s).^2);
        end
    end
    
    out = isnan(d4s{iT});
    d4s{iT}(out) = [];
    d{iT}(out,:) = [];
    p{iT}(out,:) = [];
    d4s_err{iT}(out) = [];
    spkTotalDelay{iT,1}(out,:) = [];
    spkTotalDelay{iT,2}(out,:) = [];
end

d_delay = cellfun(@(x) x(:,time>=2500 & time<=5500),d,'UniformOutput',false);
p_delay = cellfun(@(x) x(:,time>=2500 & time<=5500),p,'UniformOutput',false);
nC = cellfun(@length,d4s);
siglength = cellfun(@(x) sum(x<0.05,2),p_delay,'UniformOutput',false);



%% figures
typeList = {'IT';'PT'};
clr = {[0.745 0.102 0.125];[0.055 0.451 0.725]};
ct = cbrewer('seq','Greys',20);
ct = ct([1,4:20],:);

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 3.5]);
for iT = 1:2
    [minp,minpidx] = min(p_delay{iT},[],2);
    minp(siglength{iT}==0) = 1;
    [~,sortIdx] = sortrows([siglength{iT}>0,minpidx],{'descend','ascend'});
    
    temp = spkTotalDelay(iT,:);
    [~,maxfridx] = max(squeeze(nanmean(cat(3,temp{:}),3)),[],2);
    minpidx(minp>0.1) = 1;
    [~,sortIdx] = sortrows([siglength{iT}>0,minpidx,maxfridx],{'descend','ascend','ascend'});
    
    h(1) = axes('Position',axpt(5,1,1,1,axpt(2,1,iT,1),[0.1 0.05]));
    imagesc(1,1:nC(iT),siglength{iT}(sortIdx));
    colormap(h(1),ct);
    set(h(1),'CLim',[0 15],'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
        'XTick',[],'YTick',[1 nC(iT)],'YTickLabel',[nC(iT) 1]);
    if iT==1
        ylabel('Neurons','FontSize',8);
    end
    
    h(2) = axes('Position',axpt(5,1,2:5,1,axpt(2,1,iT,1),[0.1 0.05]));
    imagesc(time(time>=2500&time<=5500),1:nC(iT),d_delay{iT}(sortIdx,:));
    colormap(h(2),'jet')
    set(h(2),'CLim',[-1 1],'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
        'XTick',2000:1000:6000,'XTickLabel',0:1:4,'YTick',[1 nC(iT)],'YTickLabel',[]);
    xlabel('Time from delay onset (s)','FontSize',8);
    title(typeList{iT},'FontSize',9,'Color',clr{iT},'FontWeight','Bold');
end
print(fHandle,'-depsc','-painters','dprime_heatmap.ai');

%%

binRange = 0:1:sum(time>=2500&time<=5500);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 3.5]);
hold on;
for iT = 1:2
    binCounts(iT,:) = histc(siglength{iT},binRange);
    plot(binRange,cumsum(binCounts(iT,:))/sum(binCounts(iT,:)),'Color',clr{iT},'LineWidth',1);
    plot([0 0],[0 binCounts(iT,1)/sum(binCounts(iT,:))],'Color',clr{iT},'LineWidth',1);
end
xlim([-1 max(binRange)]);
ylim([0 1]);
set(gca,'Box','off','TickDir','out','FontSize',8,'XTick',0:10:30,'YTick',0:0.25:1,...
    'YTickLabel',0:25:100);
xlabel('Number of significant bins','FontSize',8);
ylabel('FON (%)','FontSize',8);
print(fHandle,'-depsc','-painters','significant_bin_num_cumsum.ai');

%%
% avedprime = cell(2,1);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 3.5]);
hold on;
for iT = 1:2
    bar(iT+0.3,nanmean(abs(d4s{iT}(siglength{iT}>0))),0.7,'EdgeColor',clr{iT},'FaceColor',clr{iT},'FaceAlpha',0.5);
    scatter(iT+rand(1,sum(siglength{iT}>0))*0.6,abs(d4s{iT}(siglength{iT}>0)),4,'o',...
        'MarkerFaceColor',clr{iT},'MarkerEdgeColor',clr{iT},'LineWidth',0.35);
    errorbar(iT+0.3,nanmean(abs(d4s{iT}(siglength{iT}>0))),nanstd(abs(d4s{iT}(siglength{iT}>0)))/...
        sqrt(sum(siglength{iT}>0)),'Color',clr{iT},'LineWidth',0.5);
end
xlim([0.5 3.1]);
ylim([0 0.6]);
set(gca,'Box','off','TickDir','out','FontSize',8,'XTick',[1.3 2.3],'YTick',0:0.2:0.6,...
    'XTickLabel',{'IT';'PT'},'XTickLabelRotation',45);
ylabel('Absolute selectivity index','FontSize',8);
pranksum = ranksum(abs(d4s{1}(siglength{1}>0)),abs(d4s{2}(siglength{2}>0)));
[~,pttest] = ttest2(abs(d4s{1}(siglength{1}>0)),abs(d4s{2}(siglength{2}>0)));
print(fHandle,'-depsc','-painters','D:\heejeong\OneDrive\git\working_memory\Figure\Nat_comm\population\abs_target_selec.ai');

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 3.5]);
hold on;
for iT = 1:2
    bar(iT+0.3,nanmean(d4s{iT}(siglength{iT}>0)),0.7,'EdgeColor',clr{iT},'FaceColor',clr{iT},'FaceAlpha',0.5);
    scatter(iT+rand(1,sum(siglength{iT}>0))*0.6,d4s{iT}(siglength{iT}>0),4,'o',...
        'MarkerFaceColor',clr{iT},'MarkerEdgeColor',clr{iT},'LineWidth',0.35);
    errorbar(iT+0.3,nanmean(d4s{iT}(siglength{iT}>0)),nanstd(d4s{iT}(siglength{iT}>0))/...
        sqrt(sum(siglength{iT}>0)),'Color',clr{iT},'LineWidth',0.5);
    psignrank(iT) = signrank(d4s{iT}(siglength{iT}>0));
    [~,pttest(iT)] = ttest(d4s{iT}(siglength{iT}>0));
end
xlim([0.5 3.1]);
ylim([-0.6 0.6]);
set(gca,'Box','off','TickDir','out','FontSize',8,'XTick',[1.3 2.3],'YTick',-0.6:0.6:0.6,...
    'XTickLabel',{'IT';'PT'},'XTickLabelRotation',45);
ylabel('Selectivity index','FontSize',8);
print(fHandle,'-depsc','-painters','target_selec.ai');


