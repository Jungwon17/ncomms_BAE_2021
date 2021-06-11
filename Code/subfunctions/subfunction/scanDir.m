function dList = scanDir(startingDir)
dList = {};
dirList = dir(startingDir);
nL = length(dirList);
if nL < 3; dList = {}; return; end;

for iL = 3:nL
    if dirList(iL).isdir == 1
        subDir = [startingDir, dirList(iL).name, '\'];
        dList = [dList; {subDir}];
        dList = [dList; scanDir(subDir)];
    end
end