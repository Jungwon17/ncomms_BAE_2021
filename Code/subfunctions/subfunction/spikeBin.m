function [time, spk] = spikeBin(spikeTime,timeWindow,binWindow,binStep)
%spikeBin Reads raster data and count spike number at trial during each time bin
%
% Syntax:  [binTime,binSpike,binIndex] = spikeBin(matFile,timeWindow,timeStep)
%   spikeTime: cell array. each cell means trial. encoded in ms.
%   timeWindow: 2 elements row vector (eq. [-1000 7000] in ms)
%   binWindow: binning window in ms (eq. 500 ms)
%   binStep: binning step in ms (eq. 100 ms)
%
% Author: Dohoung Kim
% Nov. 2015
narginchk(4,4);

bin = [timeWindow(1):binStep:(timeWindow(2)-binWindow); ...
    (timeWindow(1)+binWindow):binStep:timeWindow(2)]';
time = (timeWindow(1)+binWindow/2):binStep:(timeWindow(2)-binWindow/2);
nBin = size(bin,1);
nTrial = size(spikeTime,1);
nCell = length(spikeTime);

for iCell = 1:nCell
    if isempty(spikeTime{iCell}); 
        spikeTime{iCell}=NaN;
    end
end

spk = zeros(nTrial,nBin);
for iBin = 1:nBin
    spkTemp = cellfun(@(x) histc(x,bin(iBin,:)), spikeTime,'UniformOutput',false);
    spk(:,iBin) = cellfun(@(x) x(1), spkTemp);
end