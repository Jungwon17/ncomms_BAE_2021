function [vtTime, vtPosition, vtList] = vtLoad(vtFile)
% vtLoad searches vt files and load time and position data
%
%   vtTime: cell array of timestamp
%   vtPosition: cell array of position
%   Units in milisecond
%   vtList: shows list of VT1.nvt files
%
% Author: Junyeop Lee (cited Dohoung Kim's code)
% 

if nargin == 0;
    vtList = FindFiles('VT*.nvt','CheckSubdirs',1);
else
    if ~iscell(vtFile)
        disp('Input argument is wrong. It should be cell array.');
        return;
    elseif isempty(vtFile)
        vtList = FindFiles('VT*.nvt','CheckSubdirs',1);
    else
        nFolder = length(vtFile);
        vtList = cell(0,1);
        for iFolder = 1:nFolder
            if exist(vtFile{iFolder},'dir')
                vtList = [vtList; FindFiles('VT*.nvt','StartingDirectory',fileparts(vtFile{iFolder}),'CheckSubdirs',1)];
            elseif exist(vtFile{iFolder},'file') == 2
                [filePath, fileName, ext] = fileparts(vtFile{iFolder});
                if strcmp(ext,'.nvt')
                    vtList = [vtList; vtFile{iFolder}];
                end
            end
        end
    end
end
if isempty(vtList)
    disp('vt file does not exist!');
    [timestamp, position, vtList] = deal([]);
    return;
end
vtList = unique(vtList);
nVT = length(vtList);

vtTime = cell(nVT,1);
vtPosition = cell(nVT,1);

for iVT = 1:nVT
    [timestamp, position, ~, ~] = nvt2mat(vtList{iVT});
    timestamp = timestamp/1000; % Unit: msec
    zero = find(diff(timestamp)<0); % Delete error timestamp (timestamp should always increase)
    timestamp(zero-1) = [];
    vtTime{iVT} = timestamp;
    position(zero-1) = [];
    vtPosition{iVT} = position;
end




                    