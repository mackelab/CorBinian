function [a,b] = fitBetabinML(X, n, N, fitoptions, x0)

if nargin < 4
  fitoptions = [];
end
  
if nargin < 5
 x0 = [log(1); 
       log(1)];
end

if nargin < 3 || isempty(N)
 N = max(size(X));
end

if max(size(X)) ~= N && abs(sum(X) - 1) < 1e-10
  disp('assuming data to be given as probability mass function P(K)')
  disp('drawing N samples from P(K)')
  [~,X] = histc(rand(N,1), cumsum(X));
end
  
par = minFunc(@betabinML, x0, fitoptions, X, n);
a = exp(par(1)); b = exp(par(2));

end