function reg = slideReg(time, spk, predictor)
nBin = size(time,2);
nVar = size(predictor,2);

p = zeros(nVar,nBin);
src = zeros(nVar,nBin);
sse = zeros(nVar,nBin);

predStd = zeros(nVar,nBin);
for iVar = 1:nVar
    predStd(iVar,:) = std(predictor(:,iVar)) ./ std(spk,0,1);
end

for iBin = 1:nBin
    [beta,~,stats] = glmfit(predictor, spk(:,iBin));
    
    src(:,iBin) = beta(2:end) .* predStd(:,iBin);
    sse(:,iBin) = stats.se(2:end) .* predStd(:,iBin);
    p(:,iBin) = stats.p(2:end);
end

outRange = (abs(src+1.96*sse) > 10);
src(outRange) = 0;
sse(outRange) = 0;
p(outRange) = 1;

timesse = [time flip(time)];
sse = [src-1.96*sse flip(src+1.96*sse,2)];

reg = struct('time',time, 'p',p, 'src',src, 'timesse',timesse, 'sse',sse);
end