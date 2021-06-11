clearvars; clc; close all;

load('D:\heejeong\OneDrive\classical_conditioning\data\cellTable_ChR2.mat');

% [tData, tList] = tLoad;
cellList = T.cellList;
[cellPath,cellName,~] = cellfun(@(x) fileparts(x),cellList,'UniformOutput',false);
tList = cellfun(@(x) regexprep(x,'D:\\heejeong\\OneDrive\\classical_conditioning\\data',...
    'D:\\heejeong\\Documents\\Cheetah_data\\Classical_conditioning\\Interneuron_tagging'),...
    cellList,'UniformOutput',false);
tList = cellfun(@(x) regexprep(x,'.mat','.t'),tList,'UniformOutput',false);
nCell = length(cellList);

ChannelValidity = [1 1 1 1];
err = false(nCell,1);

for iCell = 1:nCell
    iCell
    [tPath,~,~] = fileparts(tList{iCell});
    cd(tPath)
    ttname = strsplit(cellName{iCell},'_');
    load([tPath,'\',ttname{1},'.clusters'],'-mat','MClust_Clusters');
    [~,list,~] = cellfun(@fileparts,FindFiles([ttname{1},'*.t'],...
        'StartingDirectory',fileparts(tList{iCell})),'UniformOutput',false);
    idx = strcmp(list,cellName{iCell});
    try
    spk_idx = FindInCluster(MClust_Clusters{idx});

    cluster_idx =  [tPath,'\',ttname{1},'.ntt'];
    [~, LRatio, ID] = ClusterSeparation_lratio(spk_idx, cluster_idx, ChannelValidity);
    
    save([cellPath{iCell},'\',cellName{iCell},'.mat'],'LRatio','ID','-append');
    catch error
        disp(['### Error ', tList{iCell}]);
        err(iCell) = true;
    end
end
disp('### analysis: L-Ratio & ID calculation completed! ###')