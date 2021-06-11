function [p l] = salt2(test, base, wn, dt)
%salt2   Stimulus-associated spike latency test.
%
%   Dohoung Kim, Jan. 2016.
%   Modified Balazs Hangya's salt code (balazs.cshl@gmail.com).

edges = 0:dt:wn;
nE = length(edges);

time = [base, test];

[nTrial nB] = size(time);

nhIsi = zeros(nE, nB);
for iB = 1:nB
    hst = histc(time(:,iB), edges);
    nhIsi(:, iB) = hst / sum(hst);
end

jsd = nan(nB, nB);
for iB = 1:nB
    D1 = nhIsi(:,iB);
    for jB = iB+1:nB
        D2 = nhIsi(:,jB);
        jsd(iB,jB) = sqrt(JSdiv(D1,D2)*2);
    end
end

[p, l] = makep(jsd,nB);

function [p_value Idiff] = makep(kld,kn)
% Calculates p value from distance matrix.

pnhk = kld(1:kn-1,1:kn-1);
nullhypkld = pnhk(~isnan(pnhk));   % nullhypothesis
testkld = median(kld(1:kn-1,kn));  % value to test
sno = length(nullhypkld(:));   % sample size for nullhyp. distribution
p_value = length(find(nullhypkld>=testkld)) / sno;
Idiff = testkld - median(nullhypkld);   % information difference between baseline and test latencies

% -------------------------------------------------------------------------
function D = JSdiv(P,Q)
%JSDIV   Jensen-Shannon divergence.
%   D = JSDIV(P,Q) calculates the Jensen-Shannon divergence of the two 
%   input distributions.

% Input argument check
error(nargchk(2,2,nargin))
if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
    error('Input arguments must be probability distributions.')
end
if ~isequal(size(P),size(Q))
    error('Input distributions must be of the same size.')
end

% JS-divergence
M = (P + Q) / 2;
D1 = KLdist(P,M);
D2 = KLdist(Q,M);
D = (D1 + D2) / 2;

% -------------------------------------------------------------------------
function D = KLdist(P,Q)
%KLDIST   Kullbach-Leibler distance.
%   D = KLDIST(P,Q) calculates the Kullbach-Leibler distance (information
%   divergence) of the two input distributions.

% Input argument check
error(nargchk(2,2,nargin))
if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
    error('Input arguments must be probability distributions.')
end
if ~isequal(size(P),size(Q))
    error('Input distributions must be of the same size.')
end

% KL-distance
P2 = P(P.*Q>0);     % restrict to the common support
Q2 = Q(P.*Q>0);
P2 = P2 / sum(P2);  % renormalize
Q2 = Q2 / sum(Q2);

D = sum(P2.*log(P2./Q2));
