function [f,g] = betabinML(par, X, n)

a = exp(par(1)); 
b = exp(par(2));

N = max(size(X));

%lognchoosek = (gammaln(n+1) - gammaln(X+1) - gammaln(n+1-X));  
f = sum(betaln(X + a, n - X + b)) - N * betaln(a, b);
g = [N * (psi(0, a+b) - psi(0,a) - psi(0,a+b+n)) + sum(psi(0,     X + a)); 
     N * (psi(0, a+b) - psi(0,b) - psi(0,a+b+n)) + sum(psi(0, n - X + b))];

g = g.* exp(par(:));
 
 f = -f;
 g = -g;

end