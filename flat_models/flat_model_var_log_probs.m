function [varo,meano]= flat_model_var_log_probs(ps,logbase)
%calculate variance of log-probabilities on a flat model given only the
%count distribution, without opening up the space of all states.

if nargin==1
    logbase=2;
end
scaler=1/log(logbase);

meano=-entropy_flat_model(ps,logbase);
N=numel(ps)-1;
ps=ps(:);
ps(ps==0)=1e-100;

lognchoosek=gammaln(N+1)-gammaln((0:N)'+1)-gammaln(N-(0:N)'+1);
lognchoosek=lognchoosek*scaler;
logps=log(ps)*scaler;

square_term= ps'* (logps-lognchoosek).^2;

varo=square_term-meano^2;
