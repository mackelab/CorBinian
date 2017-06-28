function [f,g, Ex, Ex2] = computeGradient2oMaxEnt(pars, N, Ek, Ek2)

mu     = pars(1);
%lambda = -exp(pars(2));
lambda = (pars(2));

ptilda = exp( mu * (0:N) + lambda * (0:N).^2 );
Z = sum( ptilda );
p      = ptilda / Z;

f = (mu * Ek + lambda * Ek2) - log(Z);

Ex  = sum( (0:N)    .* p );
Ex2 = sum( (0:N).^2 .* p );

g = [ (Ek  - Ex)  ;
      (Ek2 - Ex2)]; 
 
f = -f;
g = -g;

end