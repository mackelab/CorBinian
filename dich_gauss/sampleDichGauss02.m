function [s gamma Lambda] = sampleDichGauss02(mu,Sigma,nsamples,N,K,already_computed,acc)

% [s gamma Lambda] = sampleDichGauss02(mu,Sigma,N,K,nsamples,already_computed,acc)
%
% 	Draws nsamples samples from a Dichotomized Gaussian distribution with 
%   temporal correlations. 
%   If we want temporal correlations in s, we cannot sample s(t,:)  
%   independently. We sample s(t,:) conditioned on s(t-K+1:t-1,:). %  
%
%   Mean mu and covariance Sigma are assumed to have Hankel structure, i.e.
%   mu    : KN x 1 vector with mu(tau), 0,...,K-1 
%   Sigma : KN x KN matrix with NxN blocks given by Sigma(tau)
%
%   If you have computed the mean gamma and covariance Lambda before,
%   you can plug it in again and set already_computed to true.
%
%		Imporant: s in {0,1}^n, use appropriate covariance matrix and mean.
%		
% 	Usage: S = sampleDichGauss01([.4,.3]',[.24 .1;.1 .21],1000)
%			
% Code from the paper: 'Generating spike-trains with specified
% correlations', Macke et al., Neural Computation
%
%

if nargin<=5
    already_computed=0;
end

if nargin<=6
    acc=10^-8;
end

mu=mu(:);

if length(mu) ~= K*N
    error('ERROR: input size for mu does not match K*N')
end
if (size(Sigma,1) ~= K*N) || (size(Sigma,2) ~= K*N)
    error('ERROR: input size for Sigma does not match K*N x K*N')
end


if already_computed==0
  [gamma, Lambda] = findLatentGaussian01(mu,Sigma,acc);
else
  gamma = mu;
  Lambda = Sigma;
end

% check whether Lambda is admissable
[~, p] = chol(Lambda);
%keyboard
if p > 0 
  warning(['Covariance matrix of the latent Gaussian has at least one negative eigenvalue. ' ...
            'Applying Higham-Correction (see help higham).'])
  Lambda = higham(Lambda,1e-10,1e5);
end

% construct sampling matrices according to sec. 2.3
C = Lambda(1:N, N+1:end);
CBi = C / Lambda(1:(K-1)*N, 1:(K-1)*N);
P = chol(Lambda(1:N, 1:N) - CBi * C');
clear A C

% generate snythetic spike train
V  = zeros(nsamples, N);
dtau = -1:-1:-K+1;
g1 = gamma(1:N)';
gK = gamma(N+1:end)'; 
for t = K:nsamples                    % first K-1 entries currently ignored
    dmu = reshape(V(t+dtau,:)', 1, N*(K-1)) - gK;
    V(t,:) = g1 + dmu * CBi' + randn(1,N) * P;
end
clear t P dmu dtau CBi gK

s = V > 0; % synthetic spike-train with desired correlations






