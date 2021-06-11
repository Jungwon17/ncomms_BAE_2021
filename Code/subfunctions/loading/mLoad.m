function mList = mLoad(mFile)
%mLoad loads T*.mat files
%
%   mList: lists of TT*.mat files
%
%   Author: Dohoung Kim
%   Version 1.0 (2016/1/13)
switch nargin
    case 0
        mList = FindFiles('T*.mat','CheckSubdirs',0); 
    case 1 
        if ~iscell(mFile) 
            disp('Input argument is wrong. It should be cell array.');
            return;
        elseif isempty(mFile)
            mList = FindFiles('T*.mat','CheckSubdirs',1);
        else
            nFolder = length(mFile);
            mList = cell(0,1);
            for iFolder = 1:nFolder
                if exist(mFile{iFolder},'file')
                    mList = [mList;FindFiles('T*.mat','StartingDirectory',fileparts(mFile{iFolder}),'CheckSubdirs',1)];
                elseif strcmp(mFile{iFolder}(end-3:end),'.mat')
                    mList = [mList;mFile{iFolder}];
                end
            end
        end
end
if isempty(mList)
    disp('Mat file does not exist!');
    return;
end