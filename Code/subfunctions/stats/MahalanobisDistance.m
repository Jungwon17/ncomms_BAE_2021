function d = MahalanobisDistance(A,B)
[n1,k1] = size(A);
[n2,k2] = size(B);
n = n1+n2;
if k1~=k2
    disp('number of columns of A and B must be the same');
else
    xDiff = mean(A)-mean(B);
    cA = Covariance(A);
    cB = Covariance(B);
    pC = n1/n*cA+n2/n*cB;
    d = sqrt(xDiff*inv(pC)*xDiff');
end

function C = Covariance(X)
[n,~] = size(X);
Xc = X-repmat(mean(X),n,1);
C = Xc'*Xc/n;
