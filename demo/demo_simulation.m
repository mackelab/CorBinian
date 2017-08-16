%% criticality experiment script

clear all

% (1) Set up the simulation, fitting & evaluation processing chain
%%-------------------------------------------------------------------------

% Settings concerning the overall retina simulation
%--------------------------------------------------------------------------
dRF = [100, 100]; % dimension of visual field, in pixels
n = 50;  % number of RGCs to be simulated
d = 20;   % number of RGCs to be included in a computational run
nRuns = 1;     % number of computational runs
pRet.Ce=eye(dRF(1)*dRF(2)); % covariance matrix for Gaussian-induced noise correlations
pRet.magnitude = 1;      % parameters governing the nonlinearity mapping
pRet.gain      = 1;      % from linear filter responses to RGC output
pRet.offset    = -2.944; % spiking probability
mode.tiling = 'random';           % arrangement of RGC centres 
mode.RF = 'center-surround DoG';  % layout of RGC receptive fields

idxC = zeros(d, nRuns);
for i = 1:nRuns
    idxC(:,i) = randsample(n,d)'; % index of d cells chosen for a run
end


% Settings concerning input images to be fed into the model retina
%--------------------------------------------------------------------------
N = 2000; % number of image frames to be generated
Nc = 500;  % chunks of images to be generated at a time (memory...)
alpha = 2; % '0' for white, '1' for pink and '2' for brown noise 
xMag  = 1000; % magnitude of Gaussian noise to generate images from

% Settings concerning the maxEnt model fitting and evaluation
%--------------------------------------------------------------------------
modelFit = 'k_pairwise'; % maxEnt model to be fit to data
nChains = 1; % number of individual MCMC chains to evaluate fitting results
nSamplesEval = 10000; % number of Gibbs samples to extimate means E[f(X)]
burnIn       =  1000; % number of first Gibbs samples to be discarded
thinning     =     1; % distance in sequence between Gibbs samples to be
                      % stored (integers >1 thin out the MCMC chain)
                                            
fitoptions = struct;
fitoptions.optTol=1e-100; 
fitoptions.progTol=1e-100; 
fitoptions.MaxIter=3000;
fitoptions.MaxFunEvals=fitoptions.MaxIter;

augmentData   = true;

results.modelFit = modelFit; results.nSamplesEval = nSamplesEval;
results.burnIn = burnIn; results.thinning=thinning;

% Parameters concerning the specific layout of the RGC receptive fields
%--------------------------------------------------------------------------
pars.thres = 0.00001;% threshold below which the RFs will be truncated to 0
Sigma = { [15,   0 ;      % template for central part of DoG filters
            0,  15 ], ... % (peak of Mexican heat)
          [20,   0 ;      % template for outer part of DoG filters
            0,  20 ]};    % (surround of Mexican hat)
SigmaDFs = 100000*[1,1];  % degrees of freeom for DoG component covariances
ON_OFF = 2*randi(2, [n,1]) - 3; % per cell whether it's an ON or an OFF RGC
hight = [0.9, 1]; % template for hights of central and outer DoG components
hightSTDs = [0.01, 0.01]; % standard deviates for hights of DoG components 
idxON  = find(ON_OFF > 0); % quick lists for finding
idxOFF = find(ON_OFF < 0); % ON and OFF cells
pars.hight = zeros(n,2); % hights of DoG component Gaussian bumps
pars.hight(:,1) = abs(normrnd(hight(1),hightSTDs(1)^2,[n,1])) .*ON_OFF;
pars.hight(:,2) = abs(normrnd(hight(2),hightSTDs(2)^2,[n,1])) .*ON_OFF;
for i = 1:n 
  pars.Sigma{i,1} = wishrnd(Sigma{1}, SigmaDFs(1))/SigmaDFs(1);
  pars.Sigma{i,2} = wishrnd(Sigma{2}, SigmaDFs(2))/SigmaDFs(2);
end

%% (2) Generate input images and perform the retina simulation
%%-------------------------------------------------------------------------
disp('Generating RGC filters and spiking data')

[sim.W, sim.RGCcen] = genFilters(dRF,n,mode,pars); 
sim.W=sparse(sim.W);

sim.spikes = zeros(n,N);
for i = 1:floor(N/Nc)
  disp([' - chunk ', num2str(i), ' out of ', num2str(floor(N/Nc))])
  x = xMag * spatialPattern([dRF(1),dRF(2),Nc], -alpha);
  tmp = retSim(x,sim.W,pRet);
  sim.spikes(:,(i-1)*Nc+1:i*Nc) = tmp.spikes;
end
clear tmp x;

%% (3) Fit the statistical model to data and evaluate the quality of fit
%%-------------------------------------------------------------------------
disp('Fitting maxEnt model')
fxTrain = sim.spikes(idxC,:);
disp('- augmenting data for ML')
if augmentData % current method to avoid infinities in the ML fitting
 fxTrain(:,end+1:end+d) = triu(ones(d,d));
 fxTrain(:,end+1)       = zeros(d,1);
end
disp('- computing initialization for minFunc')
% for comparison: evaluate independent model, i.e. J = 0,L = 0
EX = full(mean(fxTrain,2));
idxL = d*(d+1)/2 + (1:d+1);
lambdaInd = zeros(d*(d+3)/2+1,1); 
mfxInd = zeros(size(lambdaInd));
lambdaInd(1:d) = log(EX./(1-EX));
lambdaInd(lambdaInd==Inf) =   1000; % fairly hacky solution in case
lambdaInd(lambdaInd==-Inf) = -1000; % EX(k) = 0 resp. EX(k) = 1

fitoptions.lambda0 = -lambdaInd;
results.lambdaInd = lambdaInd;
mfxInd(1:d) = EX;
tmp = histc(sum(bsxfun(@lt, EX, rand([d,N])),1),0:d+1)/N;
if ~isempty(idxL)
  mfxInd(idxL) =tmp(1:d+1);
end
clear EX tmp idxL tmp

fxTrain(fxTrain>1) = 1;
disp('- starting MPF fitting')
[lambdaMPF,~,~,~,~] = fit_maxent_mpf(fxTrain',fitoptions);
lambdaMPF = -lambdaMPF; % conventions...

disp('Generating data from fitted model for evaluation')
mfxTrain = zeros(d*(d+3)/2+1,1);
for i = 1:floor(N/Nc) 
  [ifxTrain, ~] = setup_features_maxent(fxTrain(:, (i-1)*Nc+1:i*Nc)', modelFit);
   mfxTrain = mfxTrain + mean(ifxTrain, 1)'; 
end
clear ifxTrain; 
mfxTrain = mfxTrain / floor(N/Nc);
mfxEval = maxEnt_gibbs_pair_C(nSamplesEval, burnIn, lambdaMPF, d);

% visualize results

figure('Name','maxEnt fitting with MPF', 'Position', [500, 500, 1200, 400])
subplot(131)
plot(mfxTrain(1:d), mfxEval(1:d), 'k.');
title('first moments')
xlabel('est.')
ylabel('data')
axis square

subplot(132)
fDescrJ = nchoosek(1:d,2)'; 
covsx = mfxTrain((d+1):(d*(d+1)/2)) - ...              
       (mfxTrain(fDescrJ(1, :)).* mfxTrain(fDescrJ(2, :))); 
covsy = mfxEval((d+1):(d*(d+1)/2)) - ...              
       (mfxEval(fDescrJ(1, :)).* mfxEval(fDescrJ(2, :))); 
plot(covsx, covsy, 'k.')
title('covariances')
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
