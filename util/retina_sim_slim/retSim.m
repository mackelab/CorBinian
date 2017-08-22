function [out] = retSim(x, W, pars)
% input:
%  -  x: d-by-N vector, gives a sequence of N input images
%  -  W: n-by-d matrix of linear filters for each neuron i = 1,...,n
%  - pars: paremeter collection to control output statistics
%       - .Ce: n-by-n covariance matrix for noise correlations
%       - .offset:   offset     of sigmoidal nonlinearity  
%       - .gain:      gain      of sigmoidal nonlinearity
%       - .magnitude: magnitude of sigmoidal nonlinearity


[n,d] = size(W);
[~,N] = size(x);

e = mvnrnd(zeros(1,d), pars.Ce, N)'; % generate noise
y = W * (x+e); % linear filtering for generating RGC input

%e = mvnrnd(zeros(1,n), pars.Ce, N)'; % generate noise
%y = y + e; % add noise on level of RGC output

% compute RGC output firing probability
y = pars.magnitude ./ (1 + exp(-(y+pars.offset)/pars.gain)); 
%y = 1 ./ (1 + exp(-y)); 

out.spikes = sparse((rand([n,N])<y));

end