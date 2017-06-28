function [gamma,lambda]=Rho_2_Lambda(mu,rho);
%function [gamma,lambda]=Rho_2_Lambda(mu,rho);
%

mu(2)=mu(1);
Sigma=[1,rho;rho,1];
Sigma=corr_2_cov_01(Sigma,mu);

[gamma,Lambda]=findLatentGaussian01(mu,Sigma);
gamma=gamma(1);
lambda=Lambda(2);

