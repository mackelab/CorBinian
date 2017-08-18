function [delta_k, LL, Sigma, deltaRecurse] = iterScalingAllVKs(...
                              Efx_K,Efy_K,N,VKold,delta0,fitoptions,Sigma)
% Given population spike count features E[f_k(X)] = P(K=k) and an old
% set of parameter V(K), here termed VKold, finds an update of parameters 
% delta_k that maximizes the log-likelihood (minimizes the neg log-likel.).
% The parameter update is regularised by a Gaussian prior on V(K) that
% combines ridge regression with a smoothness assumption.

idxBad = false(size(Efy_K));
% idxBad = (Efy_K == 0)
% finite sampling issue: if E_emp[f_k(X)] > 0 for 
% the empirical distribution p_emp(X), but E[f_k(Y)] for the MCMC sample 
% for any k, 0 <= k <= n, the most recent guess of lambda_true, 
% then we are in trouble.
% The equations for this would situation try to set V(K=k) to Infinity. 
% We'll catch these cases. 
idxBad(fitoptions.maxK+2:end) = true; % second type of error catching:
                                      % avoid updating V(K) where P(K)=0
                                      % in the original data. 

% Compute covariance matrix of Gaussian prior
if nargin < 7 || isempty(Sigma)
  Sigma = diag(fitoptions.sig2_l2 * ones(fitoptions.maxK,1)) + ...
          fitoptions.sig2_sm * ...
          exp(- 0.5*(ones(fitoptions.maxK,1)*(1:fitoptions.maxK) - ...
        (1:fitoptions.maxK)'*ones(1,fitoptions.maxK)).^2/fitoptions.tau^2);
  % The above matrix is n-by-n (instead of n+1-by-n+1 for K = 0,...,N) 
  % because delta_0 = delta_k(1) = 0 is fixed
  % Next, give credit to delta_0 = 0 being fixed, so we condition on this:              
  Sigma0Rest = fitoptions.sig2_sm * ...
                   exp(- 0.5*( (1:fitoptions.maxK).^2/fitoptions.tau^2 ) ); 
  Sigma = Sigma - ...
          (Sigma0Rest'*Sigma0Rest)/(fitoptions.sig2_l2+fitoptions.sig2_sm);
end
SigmaVKold = Sigma \ VKold(1+(1:fitoptions.maxK)); % precompute for speed

maxK = []; Eexpdelta_k = []; SigmaDelta = []; LL = []; dLL = [];    
delta_k = delta0((1:fitoptions.maxK)+1);
[delta_kMax,LL,~,~] = minFunc( @LLVKs, delta0((1:fitoptions.maxK)+1), ...
                               fitoptions, ...
                               Efx_K(1:fitoptions.maxK+1), ...
                               Efy_K(1:fitoptions.maxK+1), ...
                               N, ...
                               idxBad((1:fitoptions.maxK)+1), ...
                               Sigma, SigmaVKold);

delta_k = delta0; 
delta_k(1+(1:fitoptions.maxK)) = delta_kMax;

if nargout > 3
 deltaRecurse = sanityCheck(delta_k, Efx_K, Efy_K);
end

function [LL, dLL] = LLVKs(delta_k,Efx_K,Efy_K,N,idxBad,Sigma,SigmaVKold)
 % delta_k now contains delta-Terms for K = 1, ..., maxK. delta(K) for K=0
 % is again fixed to zero. delta(K) for K > maxK shall not be updated. 
 % idxBad tells which of the V(K)-Terms, K = 1,...,maxK are also not to be
 % updated. 
 % Efx_K, Efy_K still contain P(K=k)-Terms for K = 0, N !
  maxK = length(delta_k);
  Eexpdelta_k = Efy_K(1) + Efy_K(2:maxK+1)' * exp(delta_k);
  SigmaDelta = Sigma\delta_k; % maxK-by-maxK * maxK-by1 = maxK-by-1 vector
  LL = - delta_k' * Efx_K(2:maxK+1) + ...
         log(Eexpdelta_k) + ...
         delta_k' * (0.5*SigmaDelta + SigmaVKold)/N;
 dLL = - Efx_K(2:maxK+1) + ...
        (Efy_K(2:maxK+1) .* exp(delta_k)) / Eexpdelta_k + ...
        (SigmaDelta + SigmaVKold)/N;
 dLL(idxBad) = 0; % do not update these

end

function res = sanityCheck(delta_k, Efx_K, Efy_K)
  res = log(Efx_K ./ Efy_K) + log(Efy_K' * exp(delta_k));
end

end
