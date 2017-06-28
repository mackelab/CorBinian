function hc = binHistIndep(mu,states)
% hc = binHistIndep(mu)
% 	Computes expected histogram over binary states under independence assumption
%   P(X)=PROD(P(x_i))
%
% Code from the paper: 'Generating spike-trains with specified
% correlations', Macke et al., Neural Computation


% generate all possible binary statess

if nargin==1;
n = size(mu,1);
c = 0:2^n-1;
states = zeros(n,size(c,2));

for i=n:-1:1
    idx = c>=2^(i-1);
    states(i,idx)=1;
    c(idx) = c(idx) - 2^(i-1);    
end

states = flipud(states)';
keyboard
else
end

% transform to probabilities
mu = mu/2+.5;

% find relevant probabilities for independent model
pMat = (repmat(mu,1,size(states,2)).*states) + (repmat(1-mu,1,size(states,2)).* (~states));

% calculate histogram
hc = prod(pMat)';