function P = chisqNxN(data)
%   Pearson's Chi-Square statictics for NxN tables;
%   By Dohoung Kim, lapis42@gmail.com, Apr. 2015.
%   X_p^2 = sum {(f_ij - E_ij)^2 / E_ij}
%
%   Input:
%       indata: 2D NxN data with frequencies.
%           exam: data = [21 54; 32 46]
%
%           Y1    Y2
%   ----------------------------
%     X1    A     B     R1
%     X2    C     D     R2
%   ----------------------------
%           C1    C2    C3    W
%
%
narginchk(1,1);

C = sum(data,1);
R = sum(data,2);
data(:,C==0) = [];
C(C==0) = [];
data(R==0,:) = [];
R(R==0) = [];

nX = size(data,1);
nY = size(data,2);
if nX<=1 | nY<=1
    P = 1;
    return;
end

W = sum(data(:));

E = (R*C)/W;

chis = (data-E).^2 ./ E;

chi = sum(chis(:));
df = (nX-1) * (nY-1);

P = 1-chi2cdf(chi,df);