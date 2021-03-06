close all; clc;

%% Plotting
t_shift = 2;
time_plot = time_spec-t_shift;

minThre_corr = 10;
idx_plot1 = logical(prod(idxTrialCountDir(:,1:2)>=minThre_corr,2));

minThre_err = 2;
idx_plot2 = logical(prod(idxTrialCountDir(:,3:4)>=minThre_err,2));

idx_plot = idx_plot1&idx_plot2;

temp_plot = totS{1,1}-totS{1,2};
temp_plot = squeeze(nanmean(nanmean(temp_plot,2),1));
idx = ~isnan(temp_plot);

idx_plot = logical(idx_plot&idx);

for plot_i = 1:3
    if plot_i == 1 || plot_i == 2
        temp_plot = nanmean(totS{1,plot_i}(:,:,idx_plot),3);
    elseif plot_i == 3
        temp_plot = totS{1,1}(:,:,idx_plot)-totS{1,2}(:,:,idx_plot);
        temp_plot = nanmean(temp_plot,3);
    end
    
    figure('Position',[200 200 400 250]); hold on;
    imagesc(time_plot,freq_spec,temp_plot);
    
    plot([1 1]*0,[0 1]*110,'w-','LineWidth',1.5);
    plot([1 1]*2,[0 1]*110,'w-','LineWidth',1.5);
    plot([1 1]*6,[0 1]*110,'w-','LineWidth',1.5);
    
    colormap(jet); colorbar;
    
    if plot_i == 3
        caxis([-1 1]*0.5);
    else
        caxis([0 30]);
    end
    
    xlim([min(time_plot) max(time_plot)]);
    ylim([min(freq_spec) max(freq_spec)]);
    
    set(gca,'xtick',0:2:14,'ytick',0:20:100);
end

%% Correct vs Error
col = 'br';

figure('Position',[400 300 200 300]); hold on;

tot = cell(1,2);
idx_t = logical(time_plot>=2) & logical(time_plot<=6);
idx_f = logical(freq_spec>=4) & logical(freq_spec<=8);

for cor_i = 1:2
    tp = totS{1,cor_i}(idx_f,idx_t,idx_plot);
    tp = squeeze(nanmean(nanmean(tp,1),2))';
    nn = sum(~isnan(tp));

    tot{1,cor_i} = tp;
    
    h = bar(cor_i,nanmean(tp));
    h.FaceColor = 'none'; h.EdgeColor = col(cor_i); h.LineWidth = 1.5;
    
    h = errorbar(cor_i,nanmean(tp),nanstd(tp)/sqrt(nn),'k.-');
    h.Marker = 'none';
end

set(gca,'xtick',[]);
xtickangle(45);
xlim([0.2 2.8]); ylim([10 12.0])

[h,p] = ttest(tot{1,1},tot{1,2});
ylim auto
