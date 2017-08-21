
function [MoEs, Ts] = moments_of_energy(lambda, n, burnIn, nSamples, Ts)
% This serves to systematically test the length of MCMC samples needed to
% achieve stable estimates of the varience of energy.
% might want to discard the very first estimate Efy(:,1), Moe(:,1), as they
% might have been obtained without any burnin.

if nargin < 3
    burnIn= 100;
end
if nargin < 4
    nSamples = 1000;
end
if nargin < 5 || isempty(Ts)
    %Ts = 0.8:0.05:1.5;
    Ts = [0.8000, 0.8563, 0.9000, 0.9363, 0.9675, 0.9938, 1.0000, ...
          1.0175, 1.0388, 1.0575, 1.0775, 1.0988, 1.1200, 1.1413, ...
          1.1650, 1.1900, 1.2150, 1.2425, 1.2713, 1.3013, 1.3338, ...
          1.3687, 1.4063, 1.4463, 1.4900, 1.5375, 1.5900, 1.6500, ...
          1.7175, 1.7950, 1.8875, 2.0000];
end

N = length(Ts);
MoEs = zeros(2,N); % moments of energy 

for i = 1:N
   T = Ts(i);
   lambdaT = lambda/T; 
   disp([num2str(i),'/',num2str(N)])
   [~, MoEs(:,i)] = maxEnt_gibbs_pair_C(nSamples, burnIn, lambdaT, n);
end

end
