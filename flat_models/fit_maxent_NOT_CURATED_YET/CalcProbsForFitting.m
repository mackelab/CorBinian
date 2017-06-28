function [Q,Z,hofx,s]=CalcProbsForFitting(lambda,featurevals,hofx)
%keyboard
if nargin==2
    hofx=0;
end


s=lambda'*featurevals+hofx;

logZ=logsumexp(s');
Z=exp(logZ);
Q=exp(s-logZ);


%keyboard
