% Hosmer-Lemeshow goodness-of-fit test
% M. Farbood, April 29, 2012
% 
% Assumes that the last column in M is the outcome.  The previous columns
% are the predictors.  
%
function [HLstat HLp HLdf contingencyM] = HLtest(M, betas, numGroups)

    % Assumes there are 10 groups by default if a value is not specified.
    if ~exist('numGroups')
        numGroups = 10;
    end

    % Note this balances out the size of the group as much as possible.  
    % Change the values in numObservationsPerGroup vector if you want to hard code
    % boundaries.
    [N cols] = size(M);
    n = floor(N/numGroups);
    leftovers = N - (n*numGroups);
    noExtrasBoundary = numGroups-leftovers;
    numvars = cols - 1;
    numObservationsPerGroup = ones(numGroups,1) * n;
    if noExtrasBoundary ~= numGroups
        numObservationsPerGroup(noExtrasBoundary+1:end) = numObservationsPerGroup(noExtrasBoundary+1:end)+1;
    end

    % This is an approximation of the groups from the article Peng et al.,
    % 2002, "An Introduction to Logistic Regression
    % Analysis and Reporting" in the Journal of Educational Research.
    % They managed to have 90% of the groups have expected frequencies of 5 or
    % higher.  It was unclear how this was possible given the groups sizes they
    % provided and keeping data sorted by expected frequencies.  
    % I did my best approximation. Uncomment this line to use this grouping for
    % the special case for Peng et al., 2002.
    numObservationsPerGroup = [21 20 20 20 20 20 19 19 19 11];

    if (cols ~= length(betas))
        fprintf('ERROR: number of values in beta vector must equal number of predictors + 1\n');
        return;
    end

    % Get expected values for all data points
    numObservations = N;
    tmpM = M;
    % Add a colum of ones to represent the constant term
    tmpM = [ones(numObservations,1) tmpM];

    % Multiply all the predictors by their beta values and sum them,
    % e.g. rowSum = B0 + B1*X1 + B2*X2 + B3*X3 + ... (for all data points)
    betaMultiplier = ones(numObservations,1) * betas;    
    rowSum = sum(betaMultiplier .* tmpM(:,1:cols),2);
    % Now obtain the expected values for each data point
    % expected value = e^(B0 + B1*X1 + B2*X2 + B3*X3 + ...) / (1+e^(B0 + B1*X1 +
    % B2*X2 + B3*X3 + ...))
    allExpectedVals = exp(rowSum)./(1+exp(rowSum));
    % Sort by expected values
    [sortedM iSort] = sort(allExpectedVals);
    % Insert column of outcomes to next to the expected values
    sortedObs = tmpM(iSort,cols+1);
    sortedM = [sortedM sortedObs];

    % Caluculate contingency table for the groups
    currIndex = 1;
    contingencyM = [];
    for j=1:numGroups
        numObservations = numObservationsPerGroup(j);
        g = sortedM(currIndex:currIndex+numObservations-1,:);

        expected = sum(g(:,1));
        observed = sum(g(:,2));

        % Add values for group to the contingency table
        contingencyM = [contingencyM; j numObservations expected observed];

        currIndex = currIndex + numObservations;
    end

    % Hosmer-Lemeshow statistic = sum over all groups of (O-E)^2/(E(1-E/n))
    HLstat = sum((contingencyM(:,4) - contingencyM(:,3)).^2./(contingencyM(:,3).*(1-contingencyM(:,3)./contingencyM(:,2))));
    HLdf = numGroups - 2;
    HLp = 1 - chi2cdf(HLstat,HLdf);

    %fprintf('Hosmer-Lemeshow chisq = %g  p = %g  df = %d  N = %d\n', HLchisq, HLp, df, N);

