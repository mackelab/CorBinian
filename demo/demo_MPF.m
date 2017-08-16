% Demo script to showcase code for MPF (minimum probability flow) for 
% fitting maximum entropy models for multivariate binary data. 
%
% Instructions:
%  - Change to the code directory and run demo_MPF.m
%  - The script will take you through the functions .
%  - After each step, the script will pause. 
%  - To continue, hit any button.
%  - Read along in the demo.m file to follow what's happening.
%

% Simulation setup
%--------------------------------------------------------------------------
d=15;                  % data dimensionality
nSamplesData  = 10000; % draw from ground-truth parameters
nSamplesEval  = 10000; % draw from paramter estimates for comparison
burnIn        =  1000;

model = 'k_pairwise';

h=0.25*randn(d,1)-3.5;           % generate random bias terms h
J= 0.45*(randn(d)); J=triu(J,1); % generate interaction terms J 
lambda=hJ2lambda(h,J);           % vectorize
V = [0; linspace(3.5, -1, d)'];
lambdaTrue = [lambda;V];         % append population count terms V

%% generate training data
%--------------------------------------------------------------------------
% Initialize training data-generating MCMC chain with a sample drawn from
% a nested model (only h = lamdbdaTrue(1:d), i.e. no J, no V)
EX = exp(lambdaTrue(1:d))./(1+exp(lambdaTrue(1:d))); 
x0 = double(rand(d,1)<EX);                           

thinning     =     1; % distance in sequence between Gibbs samples to be
                      % stored (integers >1 thin out the MCMC chain)
disp('Generating training data')
xTrain = maxEnt_gibbs(nSamplesData, burnIn, thinning, lambdaTrue, x0, ...
                      model);
mfxTrain = full(mean(setup_features_maxent(xTrain', model),1))';

pause;

%% train model
%--------------------------------------------------------------------------
disp('Fitting maxEnt model')

% initialize optimization with independent model
EX = mean(xTrain,2);
lambdaInd = zeros(size(lambdaTrue)); 
lambdaInd(1:d) = log(EX./(1-EX));
lambdaInd(lambdaInd==Inf) =   1000; % catch cases EX(k)=0 resp. EX(k)=1
lambdaInd(lambdaInd==-Inf) = -1000; % 
fitoptions.lambda0 = - lambdaInd;
clear EX
                                            
fitoptions.optTol=1e-100; 
fitoptions.progTol=1e-100; 
fitoptions.MaxIter=3000;
fitoptions.MaxFunEvals=fitoptions.MaxIter;

disp('- starting MPF fitting')
[lambdaHat,~,~,~,~] = fit_maxent_mpf(xTrain',fitoptions);
lambdaHat = -lambdaHat; % conventions...

%% validate model
%--------------------------------------------------------------------------

% for small systems, we can compute P( X | lambdaTrue ) analytically
if d < 20
 [features,description,x]=setup_features_maxent(d,model);
 [~,~,Ptrue, ~]=logPMaxEnt(features,lambdaTrue);
 EX = sum(bsxfun(@times, x', Ptrue'),2);
 description(isnan(description)) = d+1;
 x1 = x; 
 x1(:,end+1) = 1; 
 EXX = sum(bsxfun(@times, (x1(:,description(1,d+1:d*(d+1)/2))...
                        .* x1(:,description(2,d+1:d*(d+1)/2)))',Ptrue'),2);
 EK = zeros(length(V),1);
 for k = 1:length(EK)
  EK(k) = sum((sum(x,2)==(k-1)) .* Ptrue);
 end
 mfxEval = [EX(:);EXX(:);EK(:)];
 clear x1 EX EXX EK description features Ptrue 

else % for large systems, sample (long) MCMC chain

 disp('Generating data from model fit')
 [mfxEval,~,~] = maxEnt_gibbs_pair_C(nSamplesEval, burnIn, ...
                                     lambdaHat, x0);
    
end % if d < 20  
                    
% visualize results

figure; 
subplot(131)
plot(mfxTrain(1:d), mfxEval(1:d), 'k.');
title('first moments')
xlabel('est.')
ylabel('data')
axis square

subplot(132)
plot(mfxTrain(d+1:end-d-1), mfxEval(d+1:end-d-1), 'k.')
title('second moments')
xlabel('est.')
ylabel('data')
axis square

subplot(133)
plot(mfxTrain(end-d:end), 'ko-')
hold on
plot(mfxEval(end-d:end), 'ro-')
hold off
xlabel('population count K')
ylabel('probability P(K)')
legend({'data', 'est.'})
title('population counts')
