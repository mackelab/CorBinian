function [mus, Rhos]=Lambda2RhoMult(gamma, Lambda)
%

mus=normcdf(gamma);

Rhos=ones(size(gamma));

for k=1:numel(gamma)
    Rhos(k)=(bivnor(-gamma(k),-gamma(k),Lambda(k))-mus(k)*mus(k))/mus(k)*(1-mus(k));
end
%keyboard
