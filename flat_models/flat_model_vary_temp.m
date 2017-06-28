function [count_distrib,Z,model]=flat_model_vary_temp(count_distrib,beta);
%fiven a model with specificed spike-count distribution, calculate
%spike-count distributions for different inverse temperatures beta;
%
%formula: for states x, P_b(x)= 1/Z_b exp(beta*lambda'*F(x));
%for spike-counts k, P_b(k)= 1/Z_b  (N choose k)^{1-beta}   P_1(k);

N=numel(count_distrib)-1;
lognchoosek=gammaln(N+1)-gammaln((0:N)+1)-gammaln(N-(0:N)+1);

logP=log(count_distrib);

count_distrib=zeros(1,numel(count_distrib));


logP_beta= beta*logP+ (1-beta)* lognchoosek;
P_beta=exp(logP_beta);
Z=sum(P_beta);
count_distrib=P_beta/Z;

if nargout==3
    model=flat_model_calc_stats(count_distrib);
end
