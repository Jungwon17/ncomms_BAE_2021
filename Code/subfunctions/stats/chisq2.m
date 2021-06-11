function [H,P] = chisq2(A,B,C,D,alpha)
%   Pearson's Chi-Square statictics for NxN tables;
%   By Dohoung Kim, lapis42@gmail.com, Apr. 2015.
%   X_p^2 = sum {(f_ij - E_ij)^2 / E_ij}
%
%   Input:
%       indata: 2D NxN data with frequencies.
%           exam: [21 54; 32 46]
%
%           Y1    Y2
%   ----------------------------
%     X1    A     B     R1
%     X2    C     D     R2
%   ----------------------------
%           C1    C2    C3    W
%
%
narginchk(3,5);

C1 = A+C;
C2 = B+D;
R1 = A+B;
R2 = C+D;
W = A+B+C+D;

if sum([A B C D]<=5)==0 | sum([A B C D])>1000
    E11 = C1.*R1./W;
    E12 = C2.*R1./W;
    E21 = C1.*R2./W;
    E22 = C2.*R2./W;
    
    chi11 = (A-E11).^2./E11;
    chi12 = (B-E12).^2./E12;
    chi21 = (C-E21).^2./E21;
    chi22 = (D-E22).^2./E22;
    
    chi = chi11+chi12+chi21+chi22;
    
    P = 1-chi2cdf(chi,1);
else % Fisher's exact test
    P = ((factorial(R1)./factorial(A)).*(factorial(R2)./factorial(B)).*...
        (factorial(C1)./factorial(C)).*(factorial(C2)./factorial(D)))./factorial(W);
end
H = P<alpha;
