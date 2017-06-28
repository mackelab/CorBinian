clear all
close all
clc 
%---------------
nSamples = 10000;
burnIn = 5000; 
d = 20;
lambda = [randn(d,1);randn(d*(d-1)/2,1)/sqrt(d);randn(d+1,1)];
x0 = d;

% Input formatting:
%--------------------------------------------------------------------------
if numel(x0) == 1
 d = x0;                                      % Generate x0 using E[X] from
 EX = exp(lambda(1:d))./(1+exp(lambda(1:d))); % a maxEnt model with only h,
 x0 = double(rand(d,1)<EX);                   % i.e. no parameters J, L
end
d = length(x0);

h = lambda(1:d);
J = lambda(d+1:d*(d+1)/2); % only uppder diag. entries of J, vectorized

L = lambda(end-d:end);     % remember, L(1) is for K=0


% Sharpen tools for the index battle ahead...
%--------------------------------------------------------------------------
m   = false(d*(d-1)/2,d); % m indexes in d-by-d matrices which d-1 entries
fm =  zeros(d-1, d);      % on the upper diagonal half share a particular
                          % index k, in booleans. fm = find(m), per column. 
for k = 1:d 
   tmp = false(d,d);   
   tmp(k,:) = true; tmp(:,k) = true;
   tmp = tmp(logical(tril(ones(d,d),-1)));
   m(:,k)   = tmp(:);
   fm(:,k) = find(m(:,k));
end
clear tmp
pairs=nchoosek(1:d,2); % needed to quickly compute the features of data x

% Start MCMC sampling
%--------------------------------------------------------------------------
xc = logical(x0); % current sample, will be continuously updated throughout

%[xSampled,E] = pwGibbsWrapper(nSamples,burnIn,d,xc,pairs,m,fm,h,J,L); % MEX

%%
tic
[xSampled,xc] = pwGibbsMaxEnt_malloc(int32(nSamples),int32(burnIn),int32(d),...
                                   xc, pairs-1, m, fm-1, h, J, L);
toc
figure; plot(xSampled); hold on;
%%
tmp = zeros(nSamples,1);
tmpxc = zeros(d,1);
tic
[xSampled,xc] = pwGibbsMaxEnt_malloc(int32(1),int32(burnIn),int32(d),...
                                   xc, pairs-1, m, fm-1, h, J, L);
for i = 2:nSamples
[tmp,xc] = pwGibbsMaxEnt_malloc(int32(1),int32(0),int32(d),...
                                xc, pairs-1, m, fm-1, h, J, L);
   xSampled = xSampled + tmp;
   tmpxc = tmpxc + xc;
end
xSampled = xSampled/nSamples;
xc = xc/nSamples;
toc
plot(xSampled,'r'); plot(xc, 'g'), hold off
