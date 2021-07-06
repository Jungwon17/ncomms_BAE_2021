clear; clc;
close all;

%% Set paths
path_data = ''; % path where data saved
paht_coh = ''; % path where results be saved

%% Data preparation
load('cellTable_v4.mat');
load('tag_v5.mat');

typeName = {'PT','IT','PC'};
D = [T.mouseNm,T.cellList, T.cellList, T.hyperLocation];

for ii = 1:size(D,1)
    temp = char(D(ii,3));
    D(ii,3) = temp(end-8:end-4);
    
    temp = char(D(ii,2));
    
    idx_dv = strfind(temp,'_');
    D(ii,2) = temp(idx_dv(3)+1:idx_dv(5)-1);
end

D_PT = D(tag.wsefr & tag.LFP & T.firingRate>0.05,:);
D_IT = D(tag.wsrxfp & tag.LFP & T.firingRate>0.05,:);
D_PC = D(tag.pc & tag.LFP & T.firingRate>0.05,:);

%% Consider FR>0.5
%%% PT
idxPT = tag.wsefr & tag.LFP & T.firingRate>0.05;
idxPT = T.firingRate(idxPT)>0.5;

%%% IT
idxIT = tag.wsrxfp & tag.LFP & T.firingRate>0.05;
idxIT = T.firingRate(idxIT)>0.5;

%%% PC
idxPC = tag.pc & tag.LFP & T.firingRate>0.05;
idxPC = T.firingRate(idxPC)>0.5;

%% Consider FR>0.5 -NS
% %%% PT
% idxPT = tag.nsefr & tag.LFP & T.firingRate>0.05;
% idxPT = T.firingRate(idxPT)>0.5;
% 
% %%% IT
% idxIT = tag.nsrxfp & tag.LFP & T.firingRate>0.05;
% idxIT = T.firingRate(idxIT)>0.5;
% 
% %%% PC
% idxPC = tag.fs & tag.LFP & T.firingRate>0.05;
% idxPC = T.firingRate(idxPC)>0.5;

%% Main
totFR = cell(2,3);
totBurst = cell(2,3);
totBurst_whole = cell(2,3);

for type_i = 1:3
    type_i

    if type_i == 1
        neuronData = D_PT;
    elseif type_i == 2
        neuronData = D_IT;
    elseif type_i == 3
        neuronData = D_PC;
    end
    
    type_folder = [paht_coh typeName{type_i} '\'];
    
    temp_folder = dir(type_folder);
    neuronN = length(temp_folder)-2;
    
    tempFR = nan(size(neuronData,1),2);
    tempBurst = nan(size(neuronData,1),2);
    tempBurst_whole = nan(size(neuronData,1),2);

    for f_i = 3:length(temp_folder)
        
        disp(['%%%' num2str(f_i-2) ' / ' num2str(neuronN)]);
        
        load([type_folder temp_folder(f_i).name]);
        
        cellNum = temp_folder(f_i).name;
        cellNum(1:2) = []; cellNum(end-3:end) = []; 
        cellNum = str2double(cellNum);
        
        tpoint = 3; % delay period
        tempFR(cellNum,:) = FRmat(:,tpoint);
        tempBurst(cellNum,:) = burstMat(:,tpoint);
        tempBurst_whole(cellNum,:) = burstMat_whole;

    end
    totFR{1,type_i} = tempFR(:,1);
    totFR{2,type_i} = tempFR(:,2);
    
    totBurst{1,type_i} = tempBurst(:,1);
    totBurst{2,type_i} = tempBurst(:,2);
    
    totBurst_whole{1,type_i} = tempBurst_whole(:,1);
    totBurst_whole{2,type_i} = tempBurst_whole(:,2);
end

%% Firing rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IT vs PT
cor_i = 1;
xpos = [2 1];

tot = cell(1,2);

col = 'cm';
figure('Position',[200 200 200 300]); hold on;
for type_i = 1:2
    if type_i == 1
        idxTrial = idxPT;
    elseif type_i == 2
        idxTrial = idxIT;
    end
    
    tp = totFR{cor_i,type_i}(idxTrial);
    nn = sum(~isnan(tp));
    
    tot{1,type_i} = tp;
    
    h = bar(xpos(type_i),nanmean(tp));
    h.FaceColor = 'w'; h.EdgeColor = col(type_i);
    h.LineWidth = 1.5;
    
    h = errorbar(xpos(type_i),nanmean(tp),nanstd(tp)/sqrt(nn),'k.');
    h.Marker = 'none'; h.LineWidth = 1.5;
end
set(gca,'xtick',[],'xtickLabel',{'PT';'IT'});
set(gca,'ytick',0:1:10);

xlim([0.2 2.8]);
ylim([0 4]);

[h,p] = ttest2(totFR{cor_i,1},totFR{cor_i,2});

%% PC: Correct vs Error
type_i = 3;

%%% index
idxTrial = idxPC;
idxPlotList = cell(1,3);

td = totCoh_cell{1,type_i}-totCoh_cell{2,type_i};
cc = squeeze(nanmean(nanmean(td,2),1));

if type_i == 3; idxPlotList{1,type_i} = ~isnan(cc) & idxTrial; end

tot =cell(1,2);

col = 'br';
figure('Position',[200 200 200 300]); hold on;
for cor_i = 1:2
    if type_i == 1
        idxTrial = idxPT;
    elseif type_i == 2
        idxTrial = idxIT;
    elseif type_i == 3
        idxTrial = idxPlotList{1,type_i};
    end
    
    tp = totFR{cor_i,type_i}(idxTrial);
    nn = sum(~isnan(tp));
    
    tot{1,cor_i} = tp;
    
    h = bar(cor_i,nanmean(tp));
    h.FaceColor = 'w'; h.EdgeColor = col(cor_i);
    h.LineWidth = 1.5;
    
    h = errorbar(cor_i,nanmean(tp),nanstd(tp)/sqrt(nn),'k.');
    h.Marker = 'none'; h.LineWidth = 1.5;
end
set(gca,'xtick',1:2,'xtickLabel',{'Correct';'Error'}); xtickangle(45);
set(gca,'ytick',0:1:10);

xlim([0.2 2.8]); ylim([0 3]);

[h,p] = ttest(tot{1,1},tot{1,2});

%% Burst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IT vs PT
cor_i = 1;
xpos = [2 1];

tot = cell(1,2);

col = 'cm';
figure('Position',[200 200 200 300]); hold on;
for type_i = 1:2
    
    if type_i == 1
        idxTrial = idxPT;
    elseif type_i == 2
        idxTrial = idxIT;
    elseif type_i == 3
    end
    
    tp = totBurst{cor_i,type_i}(idxTrial);
    nn = sum(~isnan(tp));
    
    tot{1,type_i} = tp;
    
    h = bar(xpos(type_i),nanmean(tp));
    h.FaceColor = 'w'; h.EdgeColor = col(type_i);
    h.LineWidth = 1.5;
    
    h = errorbar(xpos(type_i),nanmean(tp),nanstd(tp)/sqrt(nn),'k.');
    h.Marker = 'none'; h.LineWidth = 1.5;
end
set(gca,'xtick',[],'xtickLabel',{'PT';'IT'});
set(gca,'ytick',0:0.1:0.3);

xlim([0.2 2.8]);
ylim([0 0.3]);

[h,p] = ttest2(tot{1,1},tot{1,2});
    
%% PC: Correct vs Error
type_i = 3;

%%% index
idxTrial = idxPC;
idxPlotList = cell(1,3);

td = totCoh_cell{1,type_i}-totCoh_cell{2,type_i};
cc = squeeze(nanmean(nanmean(td,2),1));

tempPlot = cat(2,totBurst{1,type_i},totBurst{2,type_i});

td = diff(tempPlot,1,2);
idxnan = isnan(td);

if type_i == 3
    idxPlotList{1,type_i} = ~isnan(cc) & idxTrial & ~isnan(td);
end

tot = cell(1,2);

col = 'br';
figure('Position',[200 200 200 300]); hold on;
for cor_i = 1:2
    if type_i == 1
        idxTrial = idxPT;
    elseif type_i == 2
        idxTrial = idxIT;
    elseif type_i == 3
        idxTrial = idxPlotList{1,type_i};
    end
    
    tp = totBurst{cor_i,type_i}(idxTrial);
    nn = sum(~isnan(tp));
    
    tot{1,cor_i} = tp;
    
    h = bar(cor_i,nanmean(tp));
    h.FaceColor = 'w'; h.EdgeColor = col(cor_i);
    h.LineWidth = 1.5;
    
    h = errorbar(cor_i,nanmean(tp),nanstd(tp)/sqrt(nn),'k.');
    h.Marker = 'none'; h.LineWidth = 1.5;
end
set(gca,'xtick',[],'xtickLabel',{'Correct';'Error'}); xtickangle(45);
set(gca,'ytick',0:0.1:0.3);

xlim([0.2 2.8]);
ylim([0 0.3]);

[h,p] = ttest(tot{1,1},tot{1,2});