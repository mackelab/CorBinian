function [mucount,varcount]=meancorr_2_meanvar_count(mu,rho,N);
%given a population of N neurons with mean firing rate mu and pairwise
%cross correlation rho each, calculate pairwise covariance and correlation
%coefficients.

mucount=N *mu;
v= mu*(1-mu);
c=rho*v;
varcount=c*N*(N-1)+N *v; 