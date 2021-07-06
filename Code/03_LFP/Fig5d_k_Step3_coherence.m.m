close all; clc;

%% Data preparation
D = [T.mouseNm,T.cellList, T.cellList, T.hyperLocation];

for ii = 1:size(D,1)
    temp = char(D(ii,3));
    D(ii,3) = temp(end-8:end-4);
    
    temp = char(D(ii,2));
    idx_dv = strfind(temp,'_');
    D(ii,2) = temp(idx_dv(3)+1:idx_dv(5)-1);
end

%%% WS
D_PT = D(tag.wsefr & tag.LFP & T.firingRate>0.05,:);
D_IT = D(tag.wsrxfp & tag.LFP & T.firingRate>0.05,:);
D_PC = D(tag.pc & tag.LFP & T.firingRate>0.05,:);

%%% NS
% D_PT = D(tag.nsefr & tag.LFP & T.firingRate>0.05,:);
% D_IT = D(tag.nsrxfp & tag.LFP & T.firingRate>0.05,:);
% D_PC = D(tag.fs & tag.LFP & T.firingRate>0.05,:);

% idxPT = 1:size(D_PT,1);
% idxIT = 1:size(D_IT,1);
% idxPC = 1:size(D_PC,1);

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
idxPlotList = cell(1,3);

for type_i = 1:3
    
    if type_i == 1
        idxTrial = idxPT;
    elseif type_i == 2
        idxTrial = idxIT;
    elseif type_i == 3
        idxTrial = idxPC;
    end
    
    td = totCoh_cell{1,type_i}-totCoh_cell{2,type_i};
    cc = squeeze(nanmean(nanmean(td,2),1));
    if type_i == 3
        idxPlotList{1,type_i} = ~isnan(cc) & idxTrial;
        idxTrial = idxPlotList{1,type_i};
    else
        idxPlotList{1,type_i} = idxTrial;
    end
    
    for plot_i = 1:3
        figure('Position',[200 300 300 200]); hold on;
        if plot_i == 3
            tp = totCoh_cell{1,type_i}(:,:,idxTrial)-totCoh_cell{2,type_i}(:,:,idxTrial);
        else
            tp = totCoh_cell{plot_i,type_i}(:,:,idxTrial);
            
            cc = squeeze(nanmean(nanmean(tp,2),1));
            sum(~isnan(cc))
        end
        
        imagesc(time,freq,nanmean(tp,3)); colormap(jet); colorbar;
        
        timeLabel = [2 4 8 12];
        set(gca,'xtick',timeLabel,'xtickLabel',timeLabel-2);
        
        plot([1 1]*2,[0 1]*max(freq),'w-','LineWidth',2.0);
        plot([1 1]*4,[0 1]*max(freq),'w-','LineWidth',2.0);
        plot([1 1]*8,[0 1]*max(freq),'w-','LineWidth',2.0);
        
        xlim([2 10]);
        ylim([min(freq) max(freq)]);
        
        caxis([0 0.05]);
%         caxis([0 0.025]);

        if plot_i == 3
            caxis([-1 1]*0.03);
%             caxis([-1 1]*0.025);
        end
    end
end

%% Correct vs Error
t_shift = 2;
time_plot = time-t_shift;

idx_t = logical(time_plot>=2) & logical(time_plot<=6);

col = 'br';
figure('Position',[400 300 250 200]); hold on;
for type_i = 3
    if type_i == 1
        idxplot = idxPT;
    elseif type_i == 2
        idxplot = idxIT;
    elseif type_i == 3
        idxplot = idxPlotList{1,type_i};
    end
    
    tot = cell(1,2);
    for cor_i = [2 1]    
        tp = totCoh_cell{cor_i,type_i}(:,:,idxplot);
        tp = squeeze(nanmean(tp(:,idx_t,:),2))';
        nn = sum(sum(isnan(tp),2) ~= size(tp,2));
        
        tot{1,cor_i} = tp;
        
        meanP = nanmean(tp,1);
        sdP = nanstd(tp,0,1)/sqrt(nn)*1.96;
        
        x = freq; x = x'; y = meanP'; dy = sdP';
        
        h=fill([x;flipud(x)],[y-dy;flipud(y+dy)],'-','linestyle','none');
        h.FaceColor = col(cor_i);
        set(h,'facealpha',0.2);
        h = plot(x,y,'-','LineWidth',1.5);
        h.Color = col(cor_i);
    end
    xlim([min(freq) max(freq)]);
    ylim([-0.01 0.04]);
end

ylim([-0.01 0.03])

pp = nan(1,size(tp,2));
for p_i = 1:length(pp)
    [h,p] = ttest(tot{1,1}(:,p_i),tot{1,2}(:,p_i));
    pp(p_i) = p;
end
hh=pp<0.05;
hh = double(hh);

check = nanmean(tot{1,1},1)-nanmean(tot{1,2},1);
hh = hh.*(check>0);

hh(logical(hh==0)) = nan;

h = plot(freq, hh*0.025, 'ko');
h.MarkerFaceColor = h.Color; h.MarkerSize = 2;

%% PC: Correct vs Error
figure('Position',[400 300 150 300]); hold on;

tot = cell(1,2);

for type_i = 3
    if type_i == 1
        idxplot = idxPT;
    elseif type_i == 2
        idxplot = idxIT;
    elseif type_i == 3
        idxplot = idxPlotList{1,type_i};
    end
    idx_t = logical(time_plot>=2) & logical(time_plot<=6);
    idx_f = logical(freq>=4) & logical(freq<=8);
    
    for cor_i = 1:2
        tp = totCoh_cell{cor_i,type_i}(:,:,idxplot);
        tp = squeeze(nanmean(nanmean(tp(idx_f,idx_t,:),2),1))';
        nn = sum(~isnan(tp));
        
        tot{1,cor_i} = tp;

        h = bar(cor_i,nanmean(tp));
        h.FaceColor = 'none'; h.EdgeColor = col(cor_i); h.LineWidth = 1.5;
        
        h = errorbar(cor_i,nanmean(tp),nanstd(tp)/sqrt(nn),'k.-');
        h.Marker = 'none';
    end
end

set(gca,'xtick',[],'ytick',0:0.01:1);
xlim([0.2 2.8]); ylim([0 0.02]);
[h,p] = ttest(tot{1,1},tot{1,2});

ylim([0 0.02])

%% IT vs PT
t_shift = 2;
time_plot = time-t_shift;

idx_t = logical(time_plot>=2) & logical(time_plot<=6);

col = 'cm';
figure('Position',[400 300 250 200]); hold on;
tot = cell(1,2);
for type_i = 1:2
    if type_i == 1
        idxplot = idxPT;
    elseif type_i == 2
        idxplot = idxIT;
    elseif type_i == 3
        idxplot = idxPlotList{1,type_i};
    end
    
    for cor_i = 1
        tp = totCoh_cell{cor_i,type_i}(:,:,idxplot);
        tp = squeeze(nanmean(tp(:,idx_t,:),2))';

        nn = sum(sum(isnan(tp),2) ~= size(tp,2));

        meanP = nanmean(tp,1);
        sdP = nanstd(tp,0,1)/sqrt(nn)*1.96;
        
        x = freq; x = x'; y = meanP'; dy = sdP';
        
        h=fill([x;flipud(x)],[y-dy;flipud(y+dy)],'-','linestyle','none');
        h.FaceColor = col(type_i);
        set(h,'facealpha',0.2);
        h = plot(x,y,'-','LineWidth',1.5);
        h.Color = col(type_i);
    end
    tot{1,type_i} = tp;
    
    xlim([min(freq) max(freq)]);
end

ylim([-0.02 0.06])

pp = nan(1,size(tp,2));
for p_i = 1:length(pp)
    [h,p] = ttest2(tot{1,1}(:,p_i),tot{1,2}(:,p_i));
    pp(p_i) = p;
end

hh=pp<0.05;
hh = double(hh);

check = nanmean(tot{1,2},1)-nanmean(tot{1,1},1);
hh = hh.*(check>0);

hh(logical(hh==0)) = nan;

h = plot(freq, hh*0.055, 'ko');
h.MarkerFaceColor = h.Color; h.MarkerSize = 2;

%% IT vs PT
figure('Position',[400 300 150 300]); hold on;

xpos = [2 1];
tot = cell(1,2);

for type_i = 1:2
    if type_i == 1
        idxplot = idxPT;
    elseif type_i == 2
        idxplot = idxIT;
    elseif type_i == 3
        idxplot = idxPlotList{1,type_i};
    end
    
    idx_t = logical(time_plot>=2) & logical(time_plot<=6);
    idx_f = logical(freq>=4) & logical(freq<=8);
    
    tp = totCoh_cell{1,type_i}(:,:,idxplot);
    tp = squeeze(nanmean(nanmean(tp(idx_f,idx_t,:),2),1))';
    nn = sum(~isnan(tp));
    
    tot{1,type_i} = tp;
    
    h = bar(xpos(type_i),nanmean(tp));
    h.FaceColor = 'none'; h.EdgeColor = col(type_i); h.LineWidth = 1.5;
    
    h = errorbar(xpos(type_i),nanmean(tp),nanstd(tp)/sqrt(nn),'k.-');
    h.Marker = 'none';
end

set(gca,'xtick',[],'ytick',0:0.01:1);
xlim([0.2 2.8]); 
ylim([0 0.04])

[h,p] = ttest2(tot{1,1},tot{1,2});