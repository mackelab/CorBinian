function [x,Ltruncated] = sampleCovPoissonNoFit(g,L,N,nSamples)

n = length(g);
[a,b]=eig(L);
Ltruncated=(a*diag(max(diag(b),1e-5))*a');
Ltruncated=diag(sqrt(1./diag(Ltruncated)))*Ltruncated*diag(sqrt(1./diag(Ltruncated)));

%keyboard
G=chol(Ltruncated)';


x=zeros(n,nSamples);

for i=1:N
   x=x+(bsxfun(@plus,G*randn(n,nSamples),g)>0); 
end

