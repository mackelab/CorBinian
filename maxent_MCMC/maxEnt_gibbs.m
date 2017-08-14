function [xSampled] = maxEnt_gibbs(nSamples, burnIn, thinning, lambda, x0, model, mode)
% input:
% -nSamples: number of Gibbs samples to be generated
% -  burnIn: number of to-be-discarded samples from beginning of MCMC chain
% -thinning: number of samples to be discarded between each stored pair
% -  lambda: parameter vector for maxEnt model
% -      x0: EITHER a d-by-1 vector for initial member of the MCMC,
%            OR a single number specifing the dimensionality of data d
%            In the latter case, initial chain member will be drawn.
% -   model: string specifying the layout of the feature function for the
%            maxEnt model. This slim version of the code actually only
%            supports model = 'ising_count_l_0'. 
% -   mode:  string specifying what mode to operate in. If mode = 'default'
%            xSampled contains standard Bernoulli variables.  If however
%            mode = 'rb', the sampler follows 'Rao-Blackwelling' in
%            returning the Bernoulli probabilities p(x_k = 1)

% Input formatting:
%--------------------------------------------------------------------------
if numel(x0) == 1
 d = x0;                                      % Generate x0 using E[X] from
 EX = exp(lambda(1:d))./(1+exp(lambda(1:d))); % a maxEnt model with only h,
 x0 = double(rand(d,1)<EX);                   % i.e. no parameters J, L
end
d = length(x0);

if (rand>0.5)
 x0 = 1 - x0;
end
if nargin<7 || strcmp(mode,'default')
    mode = 1; % default
elseif strcmp(mode,'RB')
    mode = 0; % convention for Rao-Blackwellizing
else
    error('Unknown value for parameter mode')
end
   
if nargin>5 && ~strcmp(model, 'ising_count_l_0')
  disp('Warning: This is a slim version only supporting "ising_count_l_0"')
end

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
pairs=nchoosek(1:d,2); % needed to quickly compute the features of data x
p1 = zeros(d,1);

% Start MCMC sampling
%--------------------------------------------------------------------------
if mode % 'default'
 xSampled = false(d, nSamples);  % return Bernoulli (binary) variables
else    % Rao-Blackwellizing
 xSampled = ones(d, nSamples); % return Bernoulli probabilities
end
xc = logical(x0); % current sample, will be continuously updated throughout
idl = sum(xc)+1;  % current activity count, will also be updated throughout
for i = 1:thinning*nSamples+burnIn 
 ks = randperm(d); % one MCMC update equals one sweep through all d data- 
 for j = 1:d       % dimensions in random order
   k = ks(j);
  
  % compute p(x_k = 0)
   idl = idl - xc(k);  % current activity count, x(k) IGNORED                
   p0 = exp(L(idl)); 
  % compute p(x_k = 1)
   xc1 = xc; xc1(k) = true; 
   p1(k) = exp( h(k) ...
           + sum( J( fm(xc1(pairs(m(:,k),1))&xc1(pairs(m(:,k),2)),k) ) )...
           + L(idl+1)  ); % +1 because we now on top assume x(k) = 1
   p1(k)  = p1(k) / (p0 + p1(k)); % normalization step
 
  
  % Update chain  
   xc(k) = (rand(1) < p1(k));  % Draw new k-th entry
   idl = idl + xc(k);
 end
 % Store newest chain member
 if i>burnIn && mod(i-burnIn, thinning)==0
  if mode % convention: mode=true for 'default'. Is quicker than strcmp
   xSampled(:,(i-burnIn)/thinning) = xc;
  else    % Rao-Blackwellizing
   xSampled(:,(i-burnIn)/thinning) = p1; 
 end
end
 
end


