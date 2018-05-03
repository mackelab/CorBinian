function [lambda,logZ,logP,fitmeans,output]=fit_maxent_mpf(x, fitoptions)
% [lambda,logZ,logP,fitmeans,output]=fit_maxent_mpf(x, fitoptions)
%
% Fits k-pairwise maximum entropy models with minimum probability flow. 
% p(x) = exp( lambda' * f(x) ) / Z(lambda)
%
% Numerical optimization with Mark Schmidt's minFunc
% https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html
% 
% returns 
% lambda : parameter vector of maximum entropy model
% output : information from the numerical optimization 
% optional (uncomment below):
% logZ     : log normalizer maxEnt model with parameters lambda
% logP     : log probability for all possible data patterns
% fitmeans : expected values E[f(X)] under the maxEnt model with lambda

[n, d] = size(x); % n : number of data points, d : data dimensionality 

%if starting point is provided, use it:
if isfield(fitoptions,'lambda0')
    lambda= fitoptions.lambda0(:);
else
    %otherwise use zeros:
    lambda=zeros( d + d*(d-1)/2 + (d+1), 1);
               %  h +     J     +   V 
end

 % Parameter convention-conversion fun.
  J = diag(lambda(1:d));                   % pack contents of h into J
  J(logical(tril(ones(size(J)),-1))) = ... % pack contents of J into J
                              lambda(d+1:d*(d+1)/2); 
  J = (J + J')/2;   % MPF actually requires some matrix computations of
                  % J with the data matrix x, thus it is tedious
                  % to work with upper-triangular J and pack/unpack each 
                  % call and to instead simply work with n^2 entries for J. 
    
  V = lambda(end-d:end); % assumes d+1 entries, i.e. a feature for K = 0

  lambda = [J(:); V(:)]; % this is what the following MPF code 
                         % assumes the parameters to be structured as
   
 % Data preprocessing (data does not change over minFunc calls, do it once)
  data.x      = x';
  data.counts = sum(data.x,1);
  % Compute mask:
  data.mask   = ones(d,n);    
  
  % 1. Only work once with each pattern
  [xUnq,idxfUnq,idxUnq] = unique(x, 'rows'); xUnq = xUnq';
  countsUnq = data.counts(idxfUnq);
  % 2. Sort by activity count (narrows down which can differ by one bit)
 for k = 0:d         
  idk    = find(countsUnq == k);            % potential neighbours (by    
  idkp1  =     (countsUnq == k+1);          % virtue of having a count 
  xkn    = xUnq(:,idkp1);                   % difference of exactly one)
  fidkp1 = find(idkp1);
  % 3. Find actual neighbours by direct comparison of patterns
  for i = 1:length(idk) % for all patterns with count k ...
   idx = find(sum(bsxfun(@ne, xkn, xUnq(:,idk(i))))==1);% actual neighbours
  % 4. For each identified neightbour, check which bit mustn't be flipped
   for j = 1:length(idx) % for all of their neighouring patterns ...
    data.mask(xUnq(:,idk(i))~=xkn(:,idx(j)),idk(i)==idxUnq) = 0; 
    data.mask(xUnq(:,idk(i))~=xkn(:,idx(j)),fidkp1(idx(j))==idxUnq) = 0; 
   end
  end
 end
 % Numerical optimization step
  [lambda,f,exitflag,output] = minFunc( @K_dK_ising_PK, lambda, fitoptions, data );
  output.fs=f;
  output.exitflag=exitflag;
  
 % More parameter convention-conversion.  
  J = reshape(lambda(1:d^2),d,d);
  h = diag(J);
  J = 2*J(logical(tril(ones(size(J)),-1)));
  V = lambda(end-d:end);                     
  lambda = [ h(:); J(:); V(:) ];

 % The following code thends to bust memory for large (n,d)
  %weights = 0; % currently cannot set them otherwise with this method.  
  % fx = setup_features_maxent(x, 'ising_count_l_0');
  %[logP,logZ,~,fitmeans]= logPMaxEnt(fx,lambda,[],weights);
  logP = []; logZ = []; fitmeans = [];
end
  
  
