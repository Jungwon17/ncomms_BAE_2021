function boxplotJW(data,group,lineclr,faceclr,linewidth)
% orignated from HJ's code 2021-04-08
% edited by JW BAE 2021-04-09
groupList = unique(group);
nGroup = length(groupList);
hold on;
for iG1 = 1:nGroup
    iG=groupList(iG1);
    datag = sort(data(group==iG));
    up = datag(round(sum(group==iG)*0.75));
    low = datag(round(sum(group==iG)*0.25));
    
    outlier = datag>up+1.5*(up-low) | datag<low-1.5*(up-low);
    plot([iG iG],[up max(datag(~outlier))],'Color',lineclr(iG1,:),'LineWidth',linewidth);
    plot([iG iG],[min(datag(~outlier)) low],'Color',lineclr(iG1,:),'LineWidth',linewidth);
    rectangle('Position',[iG-0.3 low 0.6 up-low],'FaceColor',faceclr(iG1,:),'EdgeColor',lineclr(iG1,:),'LineWidth',linewidth);
    plot([iG-0.3 iG+0.3],repmat(nanmedian(datag),1,2),'Color',lineclr(iG1,:),'LineWidth',linewidth);
    plot([iG-0.2 iG+0.2 NaN iG-0.2 iG+0.2],[repmat(max(datag(~outlier)),1,2) NaN repmat(min(datag(~outlier)),1,2)],...
        'Color',lineclr(iG1,:),'LineWidth',linewidth);
end
end