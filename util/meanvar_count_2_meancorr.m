function [mu,rho,c]=meanvar_count_2_meancorr(mucount,varcount,N)
%given a population of N neurons with population count mean mucount and varianc varcount, 
%calculate pairwise mean and correlation
%coefficients.

mu=mucount/N;
v=mu*(1-mu);
c=(varcount-N*v)/N/(N-1);
rho=c/v;

%v= mu*(1-mu);
%c=rho*v;
%varcount=c*N*(N-1)+N *v; 