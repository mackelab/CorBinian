function [mu, Sigma]=Lambda_2_Sigma(gamma, Lambda)
%given latent gaussian, find mean and variance of resulting binary rvs

gamma=gamma(:);

variances=diag(Lambda);
gamma=gamma./sqrt(variances);
Lambda=cov_2_corr(Lambda);

mu=normcdf(gamma);

Sigma=ones(numel(gamma));
for k=1:numel(gamma)
    Sigma(k,k)=mu(k).*(1-mu(k));
    for kk=k+1:numel(gamma);
        Sigma(k,kk)=bivnor(-gamma(k),-gamma(kk),Lambda(k,kk))-mu(k)*mu(kk);
        Sigma(kk,k)=Sigma(k,kk);
    end
end

