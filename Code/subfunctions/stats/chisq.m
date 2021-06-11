function [P,ta] = chisq(indata);
%   Pearson's Chi-Square statictics for NxN tables;
%   By Dohoung Kim, lapis42@gmail.com, Apr. 2015.
%   X_p^2 = sum {(f_ij - E_ij)^2 / E_ij}
%
%   Input:
%       indata: 2D data matrix with 2 columns and N observation rows.
%
%           Y1    Y2    Y3
%   ----------------------------
%     X1   fq11  fq12  fq13   R1
%     X2   fq21  fq22  fq23   R2
%   ----------------------------
%           C1    C2    C3    W
%
%

if nargin <0 | nargin>2
    error('Chi-Square test: Matrix with observations must be of size NxN');
else
    j = size(indata,2);
    if j > 2
        error('Chi-Square test: column size should be 2');
    end
end

W = size(indata,1);

D = sortrows(indata,[1 2]);
[U, ~, ID] = unique(D,'rows');
fqtemp = tabulate(ID);
fq = fqtemp(:,2);
nbin = length(U);

R = zeros(nbin,1); % for X variables
C = zeros(nbin,1); % for Y variables
[strata1,~,ci] = unique(U(:,1));
nR = length(strata1);
for iR = 1:nR
    R(ci==iR) = sum(fq(ci==iR));
end
[strata2,~,cj] = unique(U(:,2));
nC = length(strata2);
for iC = 1:nC
    C(cj==iC) = sum(fq(cj==iC));
end

E = R.*C./W;
chi = sum(((fq-E).^2)./E);

P = 1-chi2cdf(chi,(nR-1)*(nC-1));
ta = table(fqtemp(:,1),U(:,1),U(:,2),fq,'VariableNames',{'ID','X','Y','f'});