function likelihood = bivbern_lik(bins, theta, lambda, LUTp11)
%BIVBERN_LIK find the likelihood of the observed data.
%Bernoulli process.
%   bins   : a Tx2 matrix containing one spike train in each column
%   theta  : support of the parameters theta
%   lambda : support of the parameter lambda
%   LUTp11      : lookup table containing the probability of observing two
%   coincident spikes, for every combination of the tree parameters theta1,
%   theta2 and lambda.

m00 = bins(1);
m10 = bins(2);
m01 = bins(3);
m11 = bins(4);

N1 = length(lambda);
N2 = length(theta);

p_fire = 1 - normcdf(theta);

% 1) compute the log-likelihood
logLik = zeros(N1, N1, N2);
% case 11
tmp = m11 * log(LUTp11);
tmp(isnan(tmp)) = 0;    % if m = 0 & p11 = 0 ==> likelihood = 1
logLik = logLik + tmp;
% case 10
LUTp10 = max(bsxfun(@minus, p_fire, LUTp11), 0);
tmp = m10 * log(LUTp10);
tmp(isnan(tmp)) = 0;
logLik = logLik + tmp;
% case 01
LUTp01 = permute(LUTp10, [2 1 3]);
tmp = m01 * log(LUTp01);
tmp(isnan(tmp)) = 0;
logLik = logLik + tmp;
% case 00
tmp = m00 * log(max(1 - LUTp10 - LUTp01 - LUTp11, 0));
tmp(isnan(tmp)) = 0;
logLik = logLik + tmp;

% 2) finally compute the likelihood
a = max(logLik(~isinf(logLik(:))));
likelihood = exp(logLik - a);

% 3) normalize the likelyhood
int1 = trapz(lambda, likelihood, 3);
int2 = trapz(theta, int1, 2);
norm = trapz(theta, int2, 1);
likelihood = likelihood / norm;

end

