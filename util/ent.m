function ent=ent(ps,logbase)
%calculate entropy given probability vector, for basis e or 2. If no
%argument is given, uses base E.
%WARNING: I changed the convention to basis 2 being the default!!

ps=ps(:);
ps(ps==0)=1e-100;
ps=ps/sum(ps);

ps(ps<exp(-700))=exp(-700);
if nargin==1
    logbase=2;
end
scaler=1/log(logbase);
logps=log(ps)*scaler;

logps(isnan(logps))=0;

ent=-sum(ps.*logps);
