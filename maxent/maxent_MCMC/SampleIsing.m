function [MC,r,p]=SampleIsing(h,J,numsamples,p)
%function [MC,r,p]=SampleIsing(h,J,numsamples,p)
%
%sample from ising model with given parameters h an J, using Gibbs Sampling
%
%co Jakob Macke, Jakob.Macke@gmail.com

mu=h*0;
Sig=J*0;

if nargin<=3
    p=struct;
end

p.outerls=0;
p.finalMC=numsamples;
p.verbose=0;

p.hinit=h;

p.Jinit=J;
p.init='writetofile';
%keyboard
[r,p]=TrainIsing(mu,Sig,p);


MC=r.MC;


