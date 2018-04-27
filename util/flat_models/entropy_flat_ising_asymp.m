function [H,h,J,ru,rd, pu, pd, m,v,]=entropy_flat_ising_asymp(mu, rho, logbase);
%
%asymptotic entropy of flat ising model with specified log-base:
% only relevant output: H entropy, h bias term and J interaction, rest is
% pretty much junk or at least not relevant in most cases
%

if nargin==2
    logbase=2;
end
scaler=1/log(logbase);

c=rho.*mu.*(1-mu);

s=sqrt(1-4*mu+4*mu.^2+4*c);

ru=.5+.5*s;
rd=.5-.5*s;
pu=(mu+ru-1)./(2*ru-1);
pd=1-pu;
H=scaler* (-log(ru)).*ru-scaler*log(1-ru).*(1-ru);

%else
    %use different values for ru, not just 1-rd
%    pu=c./(ru.^2-2*mu.*ru+mu.^2+c);
%   pd=1-pu;
%    rd=(mu-pu.*ru)./(1-pu);
%    rd(rd>1)=nan;
%    rd(rd<0)=nan;
%    H=pu.*(-log2(ru).*ru-log2(1-ru).*(1-ru))+pd.*(-log2(rd).*rd-log2(1-rd).*(1-rd));%
%
%end


m=ru.*pu+rd.*pd;
v=ru.^2.*pu+rd.^2.*pd-m.^2;    
    

J=(log(rd)-log(ru))/(rd-ru)*2;

h=-J/2;

%HC=log(ru/rd)^2/(1/ru/rd-J)*log2(exp(1))^2;

%HC=(4*(c+mu^2)-4*mu+1)/(1/ru/rd-J)*scaler^2*J^2/4;

    
%keyboard
