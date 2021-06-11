function [h,pvalue] = liddell(varargin)
%LIDDELL: Perform Liddell's exact test on matched pairs.
% Paired proportions have traditionally been compared using McNemar's test
% but an exact alternative due to Liddell is preferable. The usual reason
% for using approximate methods, when exact methods are available, is that
% the former are "quick" if sometimes "dirty". 
% There is no such justification in the present circumstances.
% Liddell FDK. Simplified exact analysis of case-referent studies; matched
% pairs; dichotomous exposure. Journal of Epidemiology and Community Health 1983;37:82-84.
%
% Syntax: 	liddell(x,alpha)
%      
%     Inputs:
%           X - 2x2 data matrix 
%           ALPHA (default 0.05) 
%     Outputs:
%           - Point estimate of relative risk (R')
%           - Exact confidence interval
%           - F statistics with p-value
%           - Power
%
%   Example:
% In the following example, a researcher attempts to determine if a drug
% has an effect on a particular disease. 
%
%                      Drug
%                  +         -
%             --------------------
%         +   |   101   |   59   |
% Placebo     |-------------------           
%         -   |   121   |   33   |
%             --------------------
%                                       
%
%   x=[101 59; 121 33];
%
%   Calling on Matlab the function: 
%             liddell(x)
%
%   Answer is:
%
% Liddell's exact test
% Point estimate of relative risk (R') =  2.0508
% Exact 95% confidence interval = 1.4904 to 2.8493
% F = 2.0167 p-value (two-side)= 0.00000443
% alpha = 0.0500  Zb = 2.7566  Power (2-tails) = 0.0058
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2008) Liddell’s test: Perform Liddell's exact test on
% matched pairs
% http://www.mathworks.com/matlabcentral/fileexchange/22024


%Input error handling
args=cell(varargin);
nu=numel(args);
if isempty(nu)
    error('Warning: Matrix of data is missed...')
elseif nu>2
    error('Warning: Max two input data are required')
end
default.values = {[],0.05};
default.values(1:nu) = args;
[x alpha] = deal(default.values{:});
%check if x is a 2x2 matrix
if ~isequal(size(x),[2 2])
    error('Input matrix must be a 2x2 matrix')
end
if ~all(isfinite(x(:))) || ~all(isnumeric(x(:)))
    error('Warning: all X values must be numeric and finite')
end
if ~isequal(x(:),round(x(:)))
    error('Warning: X data matrix values must be whole numbers')
end
if nu==2 %if necessary check alpha
    if ~isscalar(alpha) || ~isnumeric(alpha) || ~isfinite(alpha) || isempty(alpha)
        error('Warning: it is required a numeric, finite and scalar ALPHA value.');
    end
    if alpha <= 0 || alpha >= 1 %check if alpha is between 0 and 1
        error('Warning: ALPHA must be comprised between 0 and 1.')
    end
end
clear args default nu

%observed subjects with only one reaction
ob=diag(fliplr(x));
if ob(1)<ob(2)
    ob=flipud(ob);
end
r=ob(1); s=ob(2);
%perform liddell's test
%R=Maximum likelihood estimate of relative risk R as shown by Mantel &
%Haenszel. Statistical aspects of the analysis of data from retrospective
%studies of disease. J Natl Cancer Inst 1959; 22:719-748.
R=r/s;
%Confidence interval of R
P=1-alpha/2;
Ru=(r+1)/s*finv(P,2*(r+1),2*s);
Rl=r/(s+1)/finv(P,2*(s+1),2*r);
%Test of the hypothesis that R=1
F=r/(s+1);
pvalue=2*(1-fcdf(F,2*(s+1),2*r));
h = pvalue<alpha;

%Compute power
Za=abs(-realsqrt(2)*erfcinv(alpha));
N=sum(x(:));
p=min(ob./N);
pp=max(ob(1)/ob(2),ob(2)/ob(1));
num=abs(realsqrt(N*p*(pp-1)^2)-realsqrt(Za^2*(pp+1)));
denom=realsqrt(pp+1-p*(pp-1)^2);
Zb=num/denom;
pwr=(1-0.5*erfc(-Zb/realsqrt(2)))*2;

%display results
% disp('Liddell''s exact test')
% fprintf('Maximum likelihood estimate of relative risk =  %0.4f\n',R)
% fprintf('Exact %i%% confidence interval = %0.4f to %0.4f\n',(1-alpha)*100,Rl,Ru)
% fprintf('F = %0.4f p-value (two-side)= %0.8f\n',F,pvalue)
% fprintf('alpha = %0.4f  Zb = %0.4f  Power (2-tails) = %0.4f\n',alpha,Zb,pwr)