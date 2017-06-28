function [lambda,rho, gamma]=FindDGcorrelation(x)
%calculate the 'DG-correlation-coefficient lambda', the 'Pearson'-correlation rho for
%binary data x, where each row of x is a datapoint, and each column is a
%neuron/measurement dimension'
%
%JHM 08/08

mu=mean(x);
Sigma=cov(x);
rho=corrcoef(x);


[gamma,lambda]=findLatentGaussian(mu,Sigma);


if size(x,2)==2
rho=rho(2);
lambda=lambda(2);
end
