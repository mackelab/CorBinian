function x=sample_discrete(x,P,Nsamples)
%sample from discrete distribution, parametrized by N-by-d matrix x, each
%of which is a data-point, corresponding to one of the entries of P, the 
%N-by-1 vector of probabilities of each point.

Pcum=cumsum(P);

Pcum(end)=max(1,Pcum(end-1))+1e-50;
psamples=rand(Nsamples,1);

%for each sample, need to find closest index in Pcum:


[~,b] = histc(psamples,[0;Pcum]);

%keyboard

x=x(b,:);
