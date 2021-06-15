function fLFP = filt_LFP(sig1, lowerLimit, upperLimit, sF)

if nargin < 4
    sF = 2000;
end
if isempty(sF)
    sF = 2000;
end

nyquistFreq = sF/2;
lowCut = lowerLimit/nyquistFreq;
highCut = upperLimit/nyquistFreq;
filterOrder = 5;
passband = [lowCut highCut];
[Bc, Ac] = butter(filterOrder, passband);
fLFP = filtfilt(Bc, Ac, sig1);
