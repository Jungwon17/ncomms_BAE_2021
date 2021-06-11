function [tData, tList] = tLoad(tFile)
%tLoad searches t files and load data
%
%   tData: cell array. Timestamp of unit data is inside the cell array.
%   Units in millisecond.
%   tList: shows list of TTx_x.t files
%
%   Author: Dohoung Kim
%   Version 1.0 (2016/1/13)
if nargin == 0
    tList = FindFiles('T*.t','CheckSubdirs',0); 
else
    if ~iscell(tFile)
        disp('Input argument is wrong. It should be cell array.');
        return;
    elseif isempty(tFile)
        tList = FindFiles('T*.t','CheckSubdirs',1);
    else
        nFolder = length(tFile);
        tList = cell(0,1);
        for iFolder = 1:nFolder
            if exist(tFile{iFolder},'dir')
                tList = [tList;FindFiles('T*.t','StartingDirectory',tFile{iFolder},'CheckSubdirs',1)];
            elseif exist(tFile{iFolder},'file') == 2
                [pathstr, name, ext] = fileparts(tFile{iFolder});
                if strcmp(ext, '.t')
                    tList = [tList;tFile{iFolder}];
                elseif strcmp(ext, '.mat')
                    tList = [tList;mat2tt(tFile{iFolder})];
                end
            end
        end
    end
end
if isempty(tList)
    disp('t file does not exist!');
    [tData, tList] = deal([]);
    return;
end
% tList = unique(tList);
tLoad = LoadSpikes(tList,'tsflag','ts','verbose',0);

nT = length(tList);
tData = cell(nT,1);
for iT = 1:nT
    tData{iT} = Data(tLoad{iT})/10;
end