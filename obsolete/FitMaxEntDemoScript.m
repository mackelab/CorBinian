%% Demo Script to illustrate the "FitMaxEnt" Code Package
% code fits maximum entropy models to multivariate binary data, using
% user-specified feature functions and means. Current implementation only works for data
% for which the sample space can be summed over exactly, and the
% implementation is NOT optimized for memory efficiency, and has only been
% used for N<=15 so far.
%
% This code needs the code-package 'minFunc' to be available in the
% matlab-path. 
%
%

%% First, set up toy-problem to illustrate code on:
clear all
dim_x=10; %simulate 10 dimensional problem
Ns= 1000; %generate 1000 data-points;

h=randn(dim_x,1); %generate random bias terms;
J=randn(dim_x); J=triu(J,1)/sqrt(10); %generate random coupling terms:
%convert h and J to parameter lambda-- use convention that first n entries
%of lambda is a copy of h, and lower-triangular entries in J are ignored,
%and upper-triangular entries are first taken along the first row, then
%second row, etc:
lambda=hJ2lambda(h,J);

        


%now, set up parameters of code: For the moment, are extremely conservative
%on convergence criteria:
%options for parameter learning: see minFunc for their meaning
fitoptions.optTol=1e-50; 
fitoptions.progTol=1e-50; 
fitoptions.display='off';
fitoptions.MaxIter=3000;
        
                
%set up second-order feature space for binary model:
[features,description,x]=SetupFeaturesMaxEnt(dim_x,2);
            
%calculate corresponding probabilities from ground-truth maximum entropy
%model:
[logPtrue,logZtrue,Ptrue, means_true]=PMaxEnt(features,lambda);

%now, generate synthetic data from this distribution:
features_sampled=SampleDiscrete(features,Ptrue,max(Ns));

%calculate their means as we will need this as intput for the fitting
%procedure:
means_sampled=mean(features_sampled,1); clear featuressampled
                    
           
%As a small hack, we truncate means to avoid zero-counts. One could be more
%sophisticated about this, but at the moment, I am not:
minmean=1e-5;
truncator=zeros(1,numel(means_sampled));
truncator(1:dim_x)=minmean;
truncator(dim_x+1:end)=minmean^2;
means_sampled_truncated=max(means_sampled, truncator);
means_sampled_truncated=min(means_sampled_truncated, 1-truncator);
 
 
 [lambda_learned,logZlearned, Plearned, means_learned,output]=FitMaxEntLinear(features,means_sampled_truncated, fitoptions);
 Plearned=exp(Plearned);
 EntropyLogE=ent((Plearned));
          
 close all
 subplot(1,2,1)
 loglog(Plearned,Ptrue,'.')
 xlabel('Learned P')
 ylabel('True P')
 eqline
 
 subplot(1,2,2)
 plot(means_learned,means_sampled,'.')
 xlabel('learned means');
 ylabel('sampled means');
 eqline
 
 
 
 
 