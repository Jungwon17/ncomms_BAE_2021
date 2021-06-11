function varargout = logrank(X,Y)
%LOGRANK
%   p = logrank(X,Y);
%   [p,time,H(t) of group X,H(t) of group Y] = logrank(X,Y);
%
%   Example: X = [0.2 0; 0.5 1; 0.8 0; 1.2 0];
%       First column: duration of observation
%       Second column: censoring (1: if censored observation)
%   Author: Dohoung Kim
%   Version 1.0 (2015/1/9)
%% Start
n1 = size(X,1);
n2 = size(Y,1);
if n1~=0 && n2~=0
    t = tabulate([X(:,1);Y(:,1)]);
    t(t(:,3)==0,:)=[];

    m2 = tabulate(Y(:,1));
    m2(m2(:,3)==0,:)=[];
    [~,m2at,~] = intersect(t(:,1),m2(:,1));

    % m1: event/death from group 1
    % m2: event/death from group 2
    % n1: numbers in risk from group 1
    % n2: numbers in risk from group 2
    % e1: expected events from group 1
    % e2: expected events from group 2
    bl = zeros(length(t(:,1)),1);
    km = table(t(:,1),t(:,2),bl,bl,bl,bl,bl,'VariableNames',{'time','m1','m2','q1','q2','n1','n2'});
    km.m1(m2at) = km.m1(m2at) - m2(:,2);
    km.m2(m2at) = m2(:,2);
    km.n1 = n1 - [0;cumsum(km.m1(1:(end-1)))];
    km.n2 = n2 - [0;cumsum(km.m2(1:(end-1)))];

    q1 = tabulate(X(X(:,2)==1,1));
    if ~isempty(q1)
        q1(q1(:,3)==0,:)=[]; 
        [~,q1at,~] = intersect(t(:,1),q1(:,1));
        km.m1(q1at) = km.m1(q1at) - q1(:,2);
        km.q1(q1at) = q1(:,2);
    end
    q2 = tabulate(Y(Y(:,2)==1,1));
    if ~isempty(q2)
        q2(q2(:,3)==0,:)=[];
        [~,q2at,~] = intersect(t(:,1),q2(:,1));
        km.m2(q2at) = km.m2(q2at) - q2(:,2);
        km.q2(q2at) = q2(:,2);
    end
    km(km.m1==0 & km.m2==0,:) = [];
    km.p1 = 1 - km.m1./km.n1;
    km.p2 = 1 - km.m2./km.n2;
    km.H1 = 1-cumprod(km.p1);
    km.H2 = 1-cumprod(km.p2);
    km.e1 = km.n1./(km.n1+km.n2).*(km.m1+km.m2);
    km.e2 = km.n2./(km.n1+km.n2).*(km.m1+km.m2);
    km.oe1 = km.m1 - km.e1;
    v = (km.m1+km.m2).*(km.n1)./(km.n1+km.n2).*(1-(km.n1)./(km.n1+km.n2)).*((km.n1+km.n2)-(km.m1+km.m2))./(km.n1+km.n2-1);
    ch = sum(km.oe1)^2 / nansum(v);
    p = 1 - chi2cdf(ch,1);
    switch nargout
        case 1
            varargout{1} = p;
        case 4
            varargout{1} = p;
            varargout{2} = km.time;
            varargout{3} = km.H1;
            varargout{4} = km.H2;
        otherwise
            warning('Unexpected output.');
    end
else
        switch nargout
        case 1
            varargout{1} = 1;
        case 4
            varargout{1} = 1;
            varargout{2} = [];
            varargout{3} = [];
            varargout{4} = [];
        otherwise
            warning('Unexpected output.');
        end
end