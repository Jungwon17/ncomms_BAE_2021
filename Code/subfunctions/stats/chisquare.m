function [p, chi2stat, proportionT] = chisquare(contingencyT)

contingencyT(end+1,:) = sum(contingencyT);
contingencyT(:,end+1) = sum(contingencyT,2); % observed data

proportionT = contingencyT/contingencyT(end,end); % pooled estimate of proportion
observed = contingencyT(1:end-1,1:end-1);
observed = observed(:)';
ni = length(contingencyT(1,:))-1;
nj = length(contingencyT(:,1))-1;
dof = ni+nj-1; % degree of freedom

for i = 1:ni
    for j = 1:nj
        expected(j,i) = contingencyT(end,end)*proportionT(j,end)*proportionT(end,i);
    end
end
expected = expected(:)'; % expected counts under H0


chi2stat = sum((observed-expected).^2./expected);
p = 1-chi2cdf(chi2stat,dof);