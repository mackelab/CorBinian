function [varo,meano]= variance_log_probs(ps,logbase)
%calculate variance of log-probabilities of a distribution p. Second
%argument is the base of the log to be used, leave free for default (base
%2)


if nargin==1
    logbase=2;
end
ps=ps(:);
ps(ps==0)=1e-100;
ps=ps/sum(ps(:));
scaler=1/log(logbase);
logps=log(ps)*scaler;

meano=(logps'*ps);

varo=(ps'* (logps-meano).^2);

