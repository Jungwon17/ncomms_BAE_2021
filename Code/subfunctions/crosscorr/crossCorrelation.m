function [C, B] = crossCorrelation(t1, t2, binSize, binWidth)
% comparison test with c code
% much slower than c code!! (13.6 sec vs 0.6 sec)

B = -binWidth:binSize:binWidth;
nB = length(B);

C = zeros(nB, 1);

nT = length(t1);
for iT = 1:nT
    C = C + histc(t2, t1(iT)+B);
end

% C = C / (binSize * nT);