function [xSampled,E,xc] = maxEnt_gibbs_pair_C(nSamples, burnIn, lambda, x0)
% input:
% -nSamples: number of Gibbs samples to be generated
% -  burnIn: number of to-be-discarded samples from beginning of MCMC chain
% -  lambda: parameter vector for maxEnt model
% -      x0: EITHER a d-by-1 vector for initial member of the MCMC,
%            OR a single number specifing the dimensionality of data d
%            In the latter case, initial chain member will be drawn.

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

[xSampled,xc,E] = pwGibbsMaxEnt_cluster(int32(nSamples),int32(burnIn),...
                                        int32(d),...
                                   double(xc), pairs-1, m, fm-1, h, J, L);
%[xSampled,xc] = pwGibbsMaxEnt_malloc(int32(nSamples),int32(burnIn), ...
%                                     int32(d),...
%                                   double(xc), pairs-1, m, fm-1, h, J, L);
% E = [];       
end


