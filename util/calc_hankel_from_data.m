function [mu, Sigma] = calc_hankel_from_data(S,K)
% [mu, Sigma] = calc_hankel_from_data(S,K)
%
%  Constructs Hankel-structured mean mu and covariance matrix Sigma 
%  of order K from TxN data matrix S. 
%
%   Mean mu and covariance Sigma are have Hankel structure, i.e.
%   mu    : KN x 1 vector with mu(tau), 0,...,K-1 
%   Sigma : KN x KN matrix with NxN blocks given by Sigma(tau) and 
%           diagonal blocks given by Sigma(0).
%   
%   tau = 0, ..., K-1 is the time-lag between S(t,:) and S(t-tau,:). 
%   K > 1 allows Sigma to capture time-lagged correlations in S.
%   
%   Note that T > K+1 has to hold. Ideally T >> K. 

if size(S,1) < K+2
    error('size(S,1) > K+1 has to hold.')
end
    
% construct constant mu(tau) = mu, tau = 0,...,K-1
mu = repmat(mean(S,1)',K,1); % assumes stationarity

% construct time-lagged Hankel covariance matrix (KN x KN)
SH = []; for i = 1:K, SH = [SH, S(i:end-K-1+i,:)]; end
Sigma = cov(SH); 