function [H, HK, HXgivK, Pcount]    =entropy_flat_dg(gamma, lambda, N,intsteps,logbase)
%entropy of flat dichotmized gaussian model with gaussian mean gamma and
%gaussian cross-correlation lambda, using N neurons, and integrating the
%latent state over intsteps points. 
%
%output: H is entropy (base 2 unless logbase is specified)
%HK is entropy of spike count
%HgivK is entropy given spike count
%Pcount is spike count distribution

if nargin==3
    intsteps=N+100;
end
if nargin=<4
    logbase=2;
end


Ns=0:N;

Pcount=flat_dg_count_distrib(gamma,lambda,Ns,N, intsteps);

if N>1000
    Pcount=flat_dg_count_distrib_asymp(gamma, lambda, Ns/N)/N;
    Pcount(1)=flat_dg_count_distrib(gamma,lambda,0,N, intsteps);
    Pcount(end)=flat_dg_count_distrib(gamma,lambda,N,N, intsteps);
end
Pcount=Pcount/sum(Pcount);


[HK,H,HXgivK, Hpatterns]=entropy_flat_model(Pcount,logbase);

