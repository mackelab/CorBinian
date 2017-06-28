function [H,pmf,HCount]=FlatMinEnt(Ns, mu, rho)

varo=rho*mu*(1-mu);
for k=1:numel(Ns)
    [pmf{k},fval,exitflag,output,lambda]=FitMinEntLargeN(Ns(k),mu, varo);
   [HCount(k),H(k),EntPatternsgivenModel, EntPatterns]=EntropyFromSpikeCount(pmf{k}); 
end

if numel(Ns)==1
    pmf=pmf{1};
end

HCount=HCount./Ns;
H=H./Ns;