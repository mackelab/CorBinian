%
%close all
clear all

N=334; %100 neurons
Nsamples=1000; %1000 samples
mu=0.024; %probability of spike per neuron
rho=.1; %pairwise correlation between neurons

betas=10.^[-.1:.002:.1];

[gamma,lambda,DG_probs,dg_model]=fit_flat_dg(mu,rho,N);
%[s]=sample_flat_model(DG_probs,Nsamples);
    
[lambda,maxent_probs,maxent_model]=fit_flat_maxent_model(dg_model.count_distrib);
[maxent_s]=sample_flat_model(maxent_probs,Nsamples);

[count_distribs,Zs,temp_models]=flat_model_trace_temp(dg_model.count_distrib,betas);

[count_distribs_fixed,Zs_fixed,temp_models_fixed]=flat_model_trace_temp(dg_model.count_distrib,betas,1);


%%
set(gcf,'DefaultAxesColorOrder',jet(round(numel(betas(1:10:end)))));
figure;

subplot(2,4,1)
semilogy([0:N],dg_model.count_distrib,'linewidth',2); 
title('Count distribution, original model')   

subplot(2,4,2)
semilogx(betas,[temp_models.mean])
hold on
semilogx(betas,[temp_models_fixed.mean],'r')
title('Mean as function of beta')   

subplot(2,4,3)
semilogx(betas,[temp_models.corr])
hold on
semilogx(betas,[temp_models_fixed.corr],'r')
title('Corr as function of beta') 

subplot(2,4,4)
semilogx(betas,[temp_models.entropy])
hold on
semilogx(betas,[temp_models_fixed.entropy],'r')
title('Entropy as function of beta') 

subplot(2,4,5)
semilogx(betas,[temp_models.var_log_probs])
hold on
semilogx(betas,[temp_models_fixed.var_log_probs],'r')
title('Var(log probs) as function of beta') 

subplot(2,4,6)
semilogy([0:N],vertcat(temp_models(1:10:end).count_distrib))
hold on
semilogy([0:N],dg_model.count_distrib,'linewidth',3); 
ylim([1e-20,1])

subplot(2,4,7)
semilogy([0:N],vertcat(temp_models_fixed(1:10:end).count_distrib))
hold on
semilogy([0:N],dg_model.count_distrib,'linewidth',3); 
ylim([1e-20,1])
