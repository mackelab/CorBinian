function [s,hist_s]=sample_flat_model(pcount,nsamples)
%function [s,hist_s]=sample_flat_model(pcount,nsamples)
%
% sample from a 'flat' model given the distribution of spike counts. 
% 
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
