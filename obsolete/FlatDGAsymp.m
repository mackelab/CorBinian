function [P,Z]=FlatDGAsymp(gamma,lambda,r,logit);


eta=invPhi(r);

%Pold=phi(sqrt((1-lambda)/lambda).*eta-gamma/sqrt(lambda))./phi(eta)*sqrt((1-lambda)/lambda);

%P2=-1/2/lambda*((1-2*lambda)*eta.^2-2*gamma*eta*sqrt(1-lambda)+gamma^2)+0.5*log((1-lambda)/lambda);


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

%for k=1:numel(r)
%    loceta=eta(k);
%    P2(k)=exp(-1/2*(loceta-gamma*sqrt(1-lambda)/(1-2*lambda))^2/lambda*(1-2*lambda));
%    P3(k)=exp(-1/2/lambda*((1-2*lambda)*loceta^2-2*gamma*loceta*sqrt(1-lambda)+gamma^2))*sqrt((1-lambda)/lambda);
%end
%Z=1/exp(gamma^2/(2-4*lambda))/sqrt((1-lambda)/lambda)

%P2=P2/Z;
%keyboard

P(isnan(P))=0;

function phi=phi(x)
phi=normpdf(x);


function Phi=Phi(x)
phi=normcdf(x);


function invPhi=invPhi(x)
invPhi=norminv(x);















