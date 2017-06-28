function [H,ps, rs]=MinEntAsympEntropy(mu, rho);
%at the moment, only works for mu<.5;
mmu=min(mu, 1-mu);

c=rho*mmu*(1-mmu);

crit= 1/2*mmu-mmu^2;

if c<crit
    rb=mmu+(1-mmu)*rho;
    pb=mmu/rb;
    pa=1-pb;
    ps=[pa,pb,0];
    rs=[0,rb,1];
    H=-mmu*(log2(rb)+(1/rb-1).*log2(1-rb));
else
    rs=[0,.5,1];
    pb=4*mmu*(1-mmu-(1-mmu)*rho);
    pc=mmu-1/2*pb;
    pa=1-pb-pc;
    ps=[pa ,pb ,pc];
    H=pb;
end

if mu>.5
    ps([1,3])=ps([3,1]);
end
