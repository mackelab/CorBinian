function [P,Z]=flat_dg_count_distrib_asymp(gamma,lambda,r,logit);
%asymptotic spike count distribution of dichotomized gaussian

eta=invPhi(r);


if abs(lambda-0.5)>1e-5
P=-0.5*(eta-gamma*sqrt(1-lambda)/(1-2*lambda)).^2/lambda*(1-2*lambda);
Z=-0.5*(gamma^2/(1-2*lambda)+log((1-lambda)/lambda));
else
P=gamma*eta*sqrt(2);
Z=gamma^2;;
end



P=P-Z;

%keyboard
if nargin==3 || logit==0;
  P=exp(P);
  Z=exp(Z);
end



P(isnan(P))=0;

function phi=phi(x)
phi=normpdf(x);


function Phi=Phi(x)
phi=normcdf(x);


function invPhi=invPhi(x)
invPhi=norminv(x);















