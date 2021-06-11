% analysis_LRatioID
% The function calculates and returns L-ratio and ID value.
%
% To run 'analysis_LRatioID', ClusterSepration.m and sub m-files used in the fuction are needed.
% The loading engine is fixed to 'Neuralynx' format.
% Original code came from Jeonghwan Shin, and final modification was done by Joonyeup Lee.
% Jun's version 1.0 (7. 26. 2017)

[tData, tList] = tLoad;
nCell = length(tList);

ChannelValidity = [1 1 1 1];

for iCell = 1:nCell
    disp(['### Analyzing ', tList{iCell}, '...']);
    [cellPath,cellName,~] = fileparts(tList{iCell});
    ttname = strsplit(cellName,'_');

    load([ttname{1},'.clusters'],'-mat','MClust_Clusters');
    spk_idx = FindInCluster(MClust_Clusters{str2double(ttname{2})});
    cluster_idx =  [cellPath,'\',ttname{1},'.ntt'];
    [~, LRatio, ID] = ClusterSeparation_lratio(spk_idx, cluster_idx, ChannelValidity);
    
    save([cellName,'.mat'],'LRatio','ID','-append');
end
disp('### analysis: L-Ratio & ID calculation completed! ###')