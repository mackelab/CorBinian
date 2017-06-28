function [lambda,out2]=hJ2lambda(h,J)
%        [lambda]     =hJ2lambda(h,J)
%        [ h, J ]     =hJ2lambda(lambda)
% convert parameters of Ising model in representation P(x)=1/Z*exp(h'*x+ 0.5*x'J
% x)) to P(x)=1/Z*exp(lambda' f(x)) where f(x) is a feature-representation
% as returned by the function "SetupFeaturesMaxEnt".
% Works only for binary x, as then x^2 = x. 

if nargin==2 && nargout==1     %  [lambda] = hJ2lambda(h, J)
    %convert h and J to lambda:
    J=J';
    lambda=[h(:);vec(J(tril(ones(numel(h)),-1)==1))];
%    keyboard
elseif nargin==1 && nargout==2 %  [h, J]   = hJ2lambda(lambda) 
    % work in opposite direction, i.e. convert lambda to h and J
    lambda=h;
    dimo=(-1+sqrt(1+8*numel(lambda)))/2;
    if abs(rem(dimo,1))>1e-10
        error('number of dimensions does not make sense')
    else
        h=lambda(1:dimo);
        J=zeros(dimo);
        J(tril(ones(dimo),-1)==1)=lambda(dimo+1:end);
        J=J';
        %Sigma=Sigma+Sigma';
        %Sigma=Sigma+diag(mu.*(1-mu));
        lambda=h;
        out2=J*2;
    end
end
