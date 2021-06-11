function P = MantelHaenTest(X,tail)
%   MantelHaenTest, performs Mantel-Haenszel's test for k strata of 2x2 tables.
%   Mantel and Haenszel proposed this asymptotic test based on the chi-2 distr. 
%   Assuming no three-way interaction (k independent strata).
%   Program by Steinar Thorvaldsen, steinar.thorvaldsen@uit.no, Dec. 2005. 
%   Ref.: DeltaProt toolbox at http://services.cbu.uib.no/software/deltaprot/
%   Last changes 22. Dec 2010.
%   Requires Matlab 7.1 or newer, and Matlab Statistics toolbox.

%   Input:
%         X:    data matrix (size 2x2xK) of the observed frequency cells,
%               (a,b,c,d) for each stratum k.
%         tail: desired test ('lt' or 'gt': one-tail; 'ne': two-tailed(default)).
%   Output:
%         P-value
%
%   Use:  P = MantelHaenTest(Observed,'ne') 
     
%   Each stratum k must be a 2x2 table design such as:
%
%                S     Non-S    Totals
%               ------------
%    Sample-1    a       b      r1=a+b
%    Sample-2    c       d      r2=c+d
%               ------------
%    Totals     c1=a+c  c2=b+d   N=c1+c2
%
%   (S=success; Non-S=failure).

%   The test has low power for small strata sizes, and one should limit the
%   use to situations where expected cell counts of at least 5 in most of
%   the cells in each stratum (Reis et al. 1999).

%   The null hypothesis specifies that the success probabilities are equal
%   from stratum to stratum.
%   In each stratum the two rows can be viewed as data from two independent
%   binomial distributions. If these two probabilities goes in one
%   direction in one strata, and in the oposite direction in other strata,
%   the test will have less power (less ability to detect a false null
%   hypothesis).

%   Please, use the following reference:
%   Thorvaldsen, S. , Flå, T. and Willassen, N.P. (2010) DeltaProt: a software toolbox 
%       for comparative genomics. BMC Bioinformatics 2010, Vol 11:573.
%       See http://www.biomedcentral.com/1471-2105/11/573

%   Other references:
%   Mantel, N. and Haenszel, W. (1959): Statistical aspects of the analysis of data 
%       from retrospective studies of disease. J. National Cancer Inst., 22:719-748.
%   Hollander, M. and Wolfe, D.A. (1999): Nonparametric Statistical Methods.
%       John Wiley & Sons.
%   Reis, I.M et al. (1999): Exact and asymptotic tests for homogenity in several
%       2x2 tables. Statistics in Medicine, 18, p. 893-906.

if nargin < 2
    tail ='ne'; %default two-tailed test
else
    switch tail %validate the tail parameter:
    case {'ne' 'gt' 'lt'}
        % these are ok
    otherwise            
        error('MantelHaenszelTest:UnknownTail', ...
          'The ''tail'' parameter value must be ''ne'', ''gt'', or ''lt''.');
    end %switch
end %if

[I J K] = size(X);
if I ~= 2 | J ~= 2 | K < 1
    error('MantelHaenszels test: Matrix with observations must be of size 2x2xK');
end
if any(any(any(X < 0)))
   error('MantelHaenTest expects counts that are nonnegative');
end

r1s = [X(1,1,:)+X(1,2,:)]; %vector of marginal sum 1st rows
r2s = [X(2,1,:)+X(2,2,:)]; %vector of marginal sum 2nd rows
c1s = [X(1,1,:)+X(2,1,:)]; %vector of marginal sum 1st columns
c2s = [X(1,2,:)+X(2,2,:)]; %vector of marginal sum 2nd columns
Ns = [X(1,1,:)+X(1,2,:)+X(2,1,:)+X(2,2,:)]; %vector of sum of tables 2x2
% The mean and variance are calculated from the hypergeometric distribution:
Eo = (r1s.*c1s)./Ns;  %vector of null means
Varo = (r1s.*c1s.*r2s.*c2s)./((Ns.^2).*(Ns-1));  %vector of null variances
MH = (sum(X(1,1,:))-sum(Eo))/sqrt(sum(Varo));  %the Mantel-Haenszel statistic

switch tail
    case {'gt' 'lt'} % 1-tail, cf Hollander & Wolfe p. 485:
      P = 1-normcdf(abs(MH));  %approximate the normal distribution for a 1-sided test  
    otherwise % 'ne' 2-tailed
      P = 1-chi2cdf(MH^2,1);  %approximate the Chi-square distribution for a 2-sided test
end 

