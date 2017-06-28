function [ms,rhos,Hs,HCs,ps,logZs]=TraceBetaCurve(p,betas,logbase);
%given spike-count distribution p, trace out curve of spike-count
%distributions ans associated quantities using inverse temperatures beta
if nargin==2
    logbase=2
end
%scaler=1/log(logbase);

n=numel(p)-1;
r=0:n;
lognchoosek=(gammaln(n+1)-gammaln(r+1)-gammaln(n-r+1)); %*scaler;


for k=1:numel(betas);
    beta=betas(k);
   
    logp=(1-beta)*lognchoosek+beta*log(p);
    logZ=logsumexp(logp');
    logZs(k)=logZ;
    Z=exp(logZ);
    pp=exp(logp-logZ);
    
    [meano,varo]=MeanVar(pp);
    ms(k)=meano/n;
    vv=ms(k)*(1-ms(k));
    rhos(k)=(varo-n*vv)/n/(n-1)/vv;
    [HCs(k),Hs(k)]=HeatCapacityFromSpikeCount(pp,logbase);
    %HCs(k)=HCs(k).*betas(k)^2;
    ps(k,:)=pp;
    %keyboard
end
