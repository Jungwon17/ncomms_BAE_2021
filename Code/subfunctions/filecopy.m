clc; clearvars; close all;

startFolder = 'E:\Cheetah_data_jw\Session_data';
folderToMove = uigetdir(startFolder, 'Choose folders to move. All mat at files in subdirectory will be moved.');

targetStartFolder = 'D:\heejeong\OneDrive\0.working_memory';
targetFolder = uigetdir(targetStartFolder, 'Choose folder where selected files will be moved.');

cellList = FindFiles('*.mat','StartingDirectory',folderToMove, 'CheckSubdirs', true);
nC = length(cellList);

prevdir = '';
for iCell = 1:nC
    [cellDir, cellNm,~] = fileparts(cellList{iCell});
    
    newCellDir = strrep(cellDir, folderToMove, targetFolder);
    
    if ~strcmp(prevdir,newCellDir)
        mkdir(newCellDir);
    end
    copyfile(cellList{iCell},newCellDir);
    prevdir = newCellDir;
end