function [P,Z]=flat_ising_count_distrib(h,J,N,logit)
%function [P,Z]=flat_ising_count_distrib(h,J,N,logit)
%
% calculate spike-count distribution of 'flat' Ising model, i.e. for which 
% P(x)= 1/Z* exp(hk+ 0.5 k^2 J) where k= sum(x)

if nargin==3
    logit=false;
end

range=[0:N];


hofx=(gammaln(N+1)-gammaln(range+1)-gammaln(N-range+1));

H=hofx+range*h+range.^2*J;


logZ=logsumexp(H');

H=H-logZ;

if logit
    P=(H);
    Z=logZ;
else
    P=exp(H);
    Z=exp(logZ);
end

