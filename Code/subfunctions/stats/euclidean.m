function d = euclidean(x,y)
% x: test set (observationXdimension)
% y: refersence set (observationXdimesnion)
u = nanmean(y);
[nObserv,~] = size(x);
d = NaN(nObserv,1);
for iO = 1:nObserv
    difference = x(iO,:)-u;
    d(iO) = sqrt(difference*difference');
end
end