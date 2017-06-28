function [s,hist_s]=SampleFlatIsingSmart(pcount,nsamples)
%sample from an ising model with constant couplings, given the distribution
%of spike counts. In fact, it will just sample from any model which has
%distribution of spike-counts "pcount", under the assumption that all
%patterns with the same number of spikes are equally likely. 
n=numel(pcount)-1;
pcount=pcount/sum(pcount);
cum_p=[0,cumsum(pcount)];

s=zeros(n,nsamples);
perms=zeros(n,nsamples);
count=zeros(nsamples,1);
ps=rand(nsamples,1);
%try
[hist_s,count_picks]=histc(ps,cum_p);
%catch
%    keyboard
%end
hist_s=hist_s(1:end-1);
count_picks=count_picks-1;
for k=1:nsamples
    perms(:,k)=randperm(n);
    s(perms(1:count_picks(k),k),k)=1;
  %  keyboard
end
s=s==1;
%keyboard
