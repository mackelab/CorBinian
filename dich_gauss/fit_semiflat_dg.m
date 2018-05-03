function [gamma,lambda]=fit_semiflat_dg(mu,var_count,acc);
%function [gamma,lambda]=fit_semiflat_dg(mu,var_count);
%
% fit a 'semiflat' dg model to population data with given mean and correlation.
% 'Semiflat' dichomotized Gaussian models have a correlation-state that is 1D 
% and all neurons have the same loading, i.e. Lambda is a correlation matrix 
% with a constant off-diagonal lambda.
%
%function is inefficient at the moment--could be made much faster 
if nargin==2
    acc=1e-5;
end

N=numel(mu);

varo=mu.*(1-mu);

gamma=norminv(mu);

%varcount is sum of variances plus sum of covariances:

sumcov_target=var_count-sum(varo(:));




%covariance between i and j is given by

%bivnor(-g(i),-g(j),cNew) - pn

c_max=sumcov(gamma,mu,1);
c_min=sumcov(gamma,mu,-1);


if (sumcov_target>(c_max + 1e-3)) || ...
        (sumcov_target<(c_min -1e-3))
    error('A flat DG model given count-variance does not exist!')
    
end


lMin = -1;
lMax = 1;

while lMax-lMin> acc
    
    lNew = (lMax + lMin) / 2;
    if sumcov_target > sumcov(gamma,mu,lNew);
        lMin = lNew;
    else
        lMax = lNew;
    end
end

lambda=lNew;


function c=sumcov(gamma,mu,lambda);

c=0;
n=numel(gamma);
for i=1:n
    for j=i+1:n
        c=c+ bivnor(-gamma(i),-gamma(j),lambda)-mu(i)*mu(j);
    end
end

c=2*c;

