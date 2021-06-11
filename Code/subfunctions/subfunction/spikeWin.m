function spikeTime = spikeWin(spikeData, eventTime, wins)
%spikeWin makes raw spikeData into eventTime aligned data
%   spikeData: raw data from MClust. Unit must be ms.
%   eventTime: each output cell will be eventTime aligned spike data. unit must be ms.
%   win: spike within windows will be included. unit must be ms.
narginchk(3, 3);

if isempty(eventTime); spikeTime = []; return; end;

nEvent = size(eventTime);
spikeTime = cell(nEvent);
win = wins;

for iEvent = 1:nEvent(1)
    for jEvent = 1:nEvent(2)
        timeIndex = [];
        if isnan(eventTime(iEvent,jEvent)); continue; end;
        
        if iscell(wins) && length(wins)==nEvent(2)
            win = wins{jEvent};
        end
        
        if iscell(wins) && length(wins)==nEvent(1)
            win = wins{iEvent};
        end
        
        [~,timeIndex] = histc(spikeData,eventTime(iEvent,jEvent)+win);
        if isempty(timeIndex); continue; end;
        spikeTime{iEvent,jEvent} = spikeData(logical(timeIndex))-eventTime(iEvent,jEvent);
    end
end


