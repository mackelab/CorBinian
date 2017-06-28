function [P,h,J,Z,Z2,Z3,alpha,gamma]=IsiSpikeCountDistribLargeN(r,n,mu,rho,logit);

[H,ru,rd, pu, pd, m,v,h,J,HC]=FlatIsiasympEntropy(mu, rho);

rminus=rd;
rplus=ru;
pminus=pd;
pplus=pu;


%keyboard

alpha=(log(pplus)-log(pminus))/(rplus-rminus);

gamma=2*(log(rminus)-log(rplus))/(rminus-rplus);
J=gamma/n;

%J=(log(rd)-log(ru))/(rd-ru)*2;

%keyboard

h=alpha/n-J/2*n;

Z2=log(n)+log(2*sqrt(2*pi/n/(1/(rminus*rplus)-J)))+(alpha*rminus)+(n*(etae(rminus)+0.5*J*n*(rminus.^2-rminus)));


Z=log(n)+log(sqrt(2*pi/n/(1/(rminus*rplus)-J)))+log(exp(alpha*rminus)+exp(alpha*rplus))+(n*(etae(rminus)+0.5*J*n*(rminus.^2-rminus)));

%keyboard

P=(alpha*r+n*(etae(r)+0.5*J*n*(r.^2-r)));
%keyboard
%P2=(lognchoosek(n,r*n)+h*r*n+0.5*r.^2.*n.^2*J);

%Z2=sum(P2);
%keyboard
Z3=logsumexp(P(:));
P=P-Z;

if nargin==5 && logit
else
    P=exp(P);
    Z=exp(Z);
    Z2=exp(Z2);
    Z3=exp(Z3);
end

%keyboard

%P=P/Z;


function etae=etae(r)
etae=-r.*log(r)-(1-r).*log(1-r);
etae(isnan(etae))=0;
