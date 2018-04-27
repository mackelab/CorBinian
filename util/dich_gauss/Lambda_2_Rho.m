function [r, rho, c]=LambdaToRho(gamma, lambda)
%convert between moments of latent gaussian (gamma, rho), and mean,
%correlation and covariance of the resulting binary

r=normcdf(gamma(1));
r(2)=normcdf(gamma(end));
%c=bivnor(gamma(1),gamma(end),lambda)-(r(1)*r(end));

c=mvncdf([0;0], -[gamma(1); gamma(end)], [1, lambda; lambda, 1])-r(1)*r(end);
rho=c/sqrt(r(1)*(1-r(1))*r(end)*(1-r(end)));


%keyboard
if numel(gamma)==1
    r=r(1);
end
