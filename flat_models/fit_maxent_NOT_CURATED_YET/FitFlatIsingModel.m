function [h,J,beta,Z,ps,exitflag,lognchoosek,output]=FitFlatIsingModel(mu,Sigma,n,acc,steps,Bx,EBx,h,J,beta)
%find the parameters of a "flat" Ising model, i.e. one with constant
%couplings, to mean mu, pairwise covariance Sigma, and n neurons. basically
%keyboard
%a wrapper for the function IterScaling
if nargin<=3 || isnan(nan)
    acc=1e-10;
end
if nargin<=4 || isnan(steps)
    steps=1000;
end
if nargin<=5 || any(isnan(Bx))
    useBx=false;
    Bx=nan;
    EBx=nan;
else
    useBx=true;
end


mucount=n*mu;
%varcount=n^2*Sigma-n*Sigma+n*mu*(1-mu);
%C=repmat(Sigma,n,n)+eye(n)*(mu*(1-mu)-Sigma);
%varcount=sum(sum(C));
varcount=n*mu*(1-mu)+n*(n-1)*Sigma;
%keyboard

if nargin>7
    h=h;
    J=J;
    hJ0=[h,J];

elseif abs(Sigma)<10e-10;
    h=log(mu)-log(1-mu);
    J=0;
    hJ0=[h,J];
  
else %if n<=100
    hJ0(1)=log((.001+mu)/(.001+1-mu));
    hJ0(2)=Sigma/mu/(1-mu)/n/2;
end
if useBx && nargin==10
    hJ0(3)=beta;
elseif useBx
    hJ0(3)=0;
end



[h,J,Z,ps,f,lognchoosek,beta,output]=IterScalingNew(mucount,varcount,n,acc,steps,'binom',hJ0,Bx,EBx);


Z=Z(end);
exitflag=abs(Z)<inf && ~isnan(Z);
