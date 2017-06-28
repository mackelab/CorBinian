%demo for running and testing 'flat' models, i.e. models which are
%completely specified by their spike-count distribution.

close all
clear all

%% Scenario 1: Sample from different flat DG models with same mean but different 
% correlation coefficients, and compare their rasters and
% count-distributions


N=100; %100 neurons
Nsamples=500; %1000 samples
mu=.05; %probability of spike of .1 per neuron
rhos=[0,.05,.1,.2,.4]; %pairwise correlations between neurons

for k=1:numel(rhos)
    [gamma(k),lambda(k),DG_probs{k},dg_model(k)]=fit_flat_dg(mu,rhos(k),N);
    [s{k}]=sample_flat_model(DG_probs{k},Nsamples);
end


h(1)=figure

for k=1:5
subplot(2,3,k)
imagesc(-s{k}), colormap gray
xlabel('time')
ylabel('neurons');
title(['Correlation ', num2str(rhos(k))]);
end

subplot(2,3,6)
for k=1:5
   semilogy([0:N],dg_model(k).count_distrib,'color',[0,0,k/5],'linewidth',2); 
   hold on
end
ylim([.000000001,.2])

%% Scenario 2: Compare these models to Ising model


for k=1:numel(rhos)
    [h(k),J(k),ising_probs{k},ising_model(k)]=fit_flat_ising_model(mu,rhos(k),N);
    [ising_s{k}]=sample_flat_model(ising_probs{k},Nsamples);
end


h(2)=figure

for k=1:5
subplot(2,3,k)
imagesc(-ising_s{k}), colormap gray
xlabel('time')
ylabel('neurons');
title(['Correlation ', num2str(rhos(k))]);
end

subplot(2,3,6)
for k=1:5
   semilogy([0:N],ising_model(k).count_distrib,'color',[0,0,k/5],'linewidth',2); 
   hold on
end
xlim([0,100]);
ylim([.000000001,.2])


%% Scenario 3: Fit a Tkacik/Okun type of 'model' to these data-- should give perfect fit!

h(3)=figure

for k=1:numel(rhos)
    [lambdas{k},maxent_probs{k},maxent_model(k)]=fit_flat_maxent_model(dg_model(k).count_distrib);
    [maxent_s{k}]=sample_flat_model(maxent_probs{k},Nsamples);
end

for k=1:5
subplot(2,3,k)
imagesc(-maxent_s{k}), colormap gray
xlabel('time')
ylabel('neurons');
title(['Correlation ', num2str(rhos(k))]);
end

subplot(2,3,6)
for k=1:5
   semilogy([0:N],maxent_model(k).count_distrib,'color',[0,0,k/5],'linewidth',2); 
   hold on
end
xlim([0,100]);
ylim([.000000001,.2])


%% Scenario 3: Calculate entropy of these models as a function of population
% size, and compare to asymptotic results

%% 
