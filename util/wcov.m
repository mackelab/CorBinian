function [c,m]=wcov(x,P)
%wheighted covariance matrix and mean, i.e. mean and covariance of x under
%the distribution P

if nargin==1
    N=size(x,1);
    P=ones(N,1)/N;
end

wx=bsxfun(@times,full(x),P);
m=sum(wx,1);

c=x'*wx-m'*m;
