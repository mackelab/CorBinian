function [xSampled] = maxEnt_gibbs_pair(nSamples, burnIn, thinning, lambda, x0, model, mode, output)
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
% - output:  string specyfying the returned output. 'samples' returns
%            a d-by-nSamples (mode='default') resp. d*(d+3)/2+1-by-nSamples
%            (mode='rb') matrix, where columns correspond to individual 
%            samples. 'means' returns a d-by-1 resp. d*(d+3)/2+1-by-1 
%            vector of sample averages for each (feature) dimension. 

% Input formatting:
%--------------------------------------------------------------------------
if numel(x0) == 1
 d = x0;                                      % Generate x0 using E[X] from
 EX = exp(lambda(1:d))./(1+exp(lambda(1:d))); % a maxEnt model with only h,
 x0 = double(rand(d,1)<EX);                   % i.e. no parameters J, L
end
d = length(x0);

%if (rand>0.5) % for debuggin purposes, can flip initial chain element
% x0 = 1 - x0; % to see the influence of starting conditions on the
%end           % results of the full MCMC chain

if nargin>5 && ~strcmp(model, 'ising_count_l_0')
  disp('Warning: This is a slim version only supporting "ising_count_l_0"')
end

if nargin<7 || strcmp(mode,'default')
    mode = 1; % default
elseif strcmp(mode,'rb')
    mode = 0; % convention for Rao-Blackwellizing
else
    error('Unknown value for parameter mode')
end

if nargin<8 || strcmp(output,'samples')
    output = 1; % default
elseif strcmp(output,'means')
    output = 0; % convention for Rao-Blackwellizing
else
    error('Unknown value for parameter output')
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
clear tmp
pairs=nchoosek(1:d,2); % needed to quickly compute the features of data x

% Start MCMC sampling
%--------------------------------------------------------------------------
if output % return all MCMC samples in a d-by-1 resp. d*(d+1)/2-by-1 matrix
 if mode % 'default'
  xSampled = false(d, nSamples);  % return Bernoulli (binary) variables
 else    % Rao-Blackwellizing   
  xSampled = zeros(d*(d+3)/2+1, nSamples); % return Bernoulli probabilities
 end
else      % only return the estimated expected values in a single vector
 if mode % 'default'
  xSampled = false(d, 1);  % return Bernoulli (binary) variables
 else    % Rao-Blackwellizing   
  xSampled = zeros(d*(d+3)/2+1, 1);      % return Bernoulli probabilities
 end   
 xTmp = xSampled;% intermediate storage for normalization (numerical issue)
end
xc = logical(x0); % current sample, will be continuously updated throughout
idl = sum(xc)+1;  % current activity count, will also be updated throughout
for i = 1:thinning*nSamples+burnIn 
 ks = randperm(d*(d-1)/2); % one MCMC update equals one sweep through all d data- 
                           % dimensions in random order
 pc = zeros(d*(d+3)/2+1,1);  % current Bernoulli probabilities
 for j = 1:length(ks)          
   k = pairs(ks(j),1); 
   l = pairs(ks(j),2);
   
   idl = idl - xc(k) - xc(l); % current activity count, x(k), xc(l) IGNORED                
  
  % compute some key values
   xc1 = xc; 
   xc1(k) = true; xc1(l) = false; 
   sk = sum( J( fm(xc1(pairs(m(:,k),1))&xc1(pairs(m(:,k),2)),k) ) );
   % m(:,k) holds which of the d-1 indices of J are involved with x(k)
   % pairs(m(:,k),1:2) gives the d-1 index pairs (i,j) of x that correspond
   %                   to the indices of J that are involved with x(k)
   % xc1(pairs(...,1))&x1(pairs(...,2)) filters out those where x(i)=x(j)=1
   % fm(xc1(...)&xc1(...)) translates the pair indices back into J indices
   
   xc1(k) = false; xc1(l) = true; 
   sl = sum( J( fm(xc1(pairs(m(:,l),1))&xc1(pairs(m(:,l),2)),l) ) );
   
   xc1(k) = true;  
   slk = sum( J( fm(xc1(pairs(m(:,k),1))&xc1(pairs(m(:,k),2)),k) ) );
   
   L1 = L(idl+1); % one of the two variables xc(k), xc(l) is = 1
   L2 = L(idl+2); % both variables are = 1
      
  % compute joint probabilities for bivariate Bernoulli   
   p = exp([h(k) +        sk  +       L1           % xc(k) = 1, xc(l) = 0
            h(k) + h(l) + sl  + slk + L2;          % xc(k) = 1, xc(l) = 1
                   h(l) +        sl + L1;          % xc(k) = 0, xc(l) = 1
                                      L(idl)]);    % xc(k) = 0, xc(l) = 0
   p = p/sum(p); % normalize
  % update chain   
   rnd = rand(1);               % essentially, use the CDF for each of the
   xc(k) = (rnd < p(1)+p(2));   % four possible outcomes for (xc(k),xc(l))
   xc(l) = (p(1) < rnd && rnd < 1-p(4));       
   idl = idl + xc(k) + xc(l);
   if ~mode
     pc(k) = pc(k) + p(1) + p(2);  % for E[x(k)]
     pc(l) = pc(l) + p(2) + p(3);  % for E[x(l)]
     pc(d+ks(j)) = p(2); % for E[x(k) * x(l)] and thus for Cov(x(k),x(l))
     pc(d*(d+1)/2+idl) = pc(d*(d+1)/2+idl)+1; % for E[sum(x)==k], i.e. V(K)
   end
 end % end for j=1:length(ks), i.e. loop over all updates within one sweep
 
 if ~mode
   pc(1:d) = pc(1:d) / (d-1);  % we just visited each x(k) d-1 times
   pc(d*(d+1)/2+(1:d+1)) = pc(d*(d+1)/2+(1:d+1))/ length(ks);
 end
 
 ii = i-burnIn;
 if output % return samples
  % Store newest chain member
  if ii>0 && mod(ii, thinning)==0
   if mode % convention: mode=true for 'default'. Is quicker than strcmp
    xSampled(:,ii/thinning) = xc;
   else    % Rao-Blackwellizing
    xSampled(:,ii/thinning) = pc; 
   end
  end
 else      % return means
  % Update online estimate of expected values
  if ii>0
   if mode % convention: mode=true for 'default'. 
    xTmp = xTmp + xc;
   else    % Rao-Blackwellizing
    xTmp = xTmp + pc; 
   end
   if ~mod(ii,1000)% normalize output every 1000 samples to avoid numerical 
    xSampled = ((ii-1000)*xSampled + xTmp)/ii; % issues with large numbers.
    xTmp = 0;
    disp(['  - ', num2str(ii),'/',num2str(nSamples),'samples'])
   end
  end
 end
      
end % end for i=1:nSamples
%xSampled = ((i-mod(i,1000)) * xSampled + xTmp)/i;

end


