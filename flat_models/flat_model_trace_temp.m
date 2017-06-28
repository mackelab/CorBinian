function [count_distribs,Zs,models,h]=flat_model_trace_temp(count_distrib,betas,fix_mean)
%fiven a model with specificed spike-count distribution, calculate
%spike-count distributions for different inverse temperatures beta;
%
%formula: for states x, P_b(x)= 1/Z_b exp(beta*lambda'*F(x));
%for spike-counts k, P_b(k)= 1/Z_b  (N choose k)^{1-beta}   P_1(k);

if nargin==2
    fix_mean=false;
end


if ~fix_mean
    for k=1:numel(betas);
        
        [count_distribs(k,:),Zs(k,:),locmodel]=flat_model_vary_temp(count_distrib,betas(k));
        locmodel.beta=betas(k);
        models(k)=locmodel;
    end
    
else
    mean_rate=calc_mean_var(count_distrib)/(numel(count_distrib)-1);
    for k=1:numel(betas);
        [h(k),logZ,count_distribs(k,:),locmodel]=flat_model_vary_temp_fix_mean(mean_rate,count_distrib,betas(k));
        locmodel.beta=betas(k);
        models(k)=locmodel;
        Zs(k)=exp(logZ);    
    end
end
