function dList = FindDir(findingDir, startingDir)

dList = scanDir(startingDir);
isWord = cellfun(@(x) ~isempty(strfind(x, findingDir)), dList);
dList = dList(isWord);
