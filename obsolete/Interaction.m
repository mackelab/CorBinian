function a=Interaction(mu,rho);

if numel(mu)==1
    mu=[mu;mu];
end

vars=mu.*(1-mu);
%first, calculate covariance:
c=rho*sqrt(prod(vars));
P11=c+prod(mu);
P10=mu(1)-P11;
P01=mu(2)-P11;
P00=1-P11-P10-P01;

a=1/4*log2(P00*P11/P01/P10);



%keyboard
