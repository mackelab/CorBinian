function [logP,logZ,P, means]=logPMaxEnt(x,lambda,logZ, weights, mode)
%calculate log probabilities, log partition function, probabilites and
%feature means under a model of the form P(x)=exp(-logZ+lambda'*x).
%If fourth argument 'weights' is added, it instead calculates P(x)=
%exp(-logZ+ lambda'*x + weights); so weights needs to be a vector of the
%length size(x,2);

if numel(x)==1
    dimo=x;
    %use standard ordering for x, assuming a second order Ising model:
    %dimo=dimo=-1+sqrt(1+4*numel(lambda));
    switch mode
        case {2,'ising'}
           [x]=setup_features_maxent(dimo,2);
        case 'ising_count'
           [x]=setup_features_maxent(dimo,'ising_count'); 
        case 'ising_count_l_0'
           [x]=setup_features_maxent(dimo,'ising_count_l_0'); 
    end
    % keyboard
end


if nargin<=3
    weights=0;
end

logQ=full(x*lambda)+weights;


if nargin<=2 || isempty(logZ)
    logZ=logsumexp(logQ);
end


logP=logQ-logZ;


if nargout>=3
    P=exp(logP);
    SumP=sum(P);
    P=P/SumP;
end

if nargout==4
    means=P(:)'*x;
end


