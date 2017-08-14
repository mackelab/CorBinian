%% criticality experiment script
function mainExperiment(d, n, nComp, nRuns, N, Nc, alpha, xMag, nChains, nSamplesEval, burnIn)
%% ------------------------------------------------------------------------
% (1) Set up the simulation, fitting & evaluation processing chain
%%-------------------------------------------------------------------------

% Settings concerning the overall retina simulation
%--------------------------------------------------------------------------
if isempty(d)
    d = [100, 100]; % dimension of visual field, in pixels
end
if isempty(n)
    n = 200;       % number of RGCs to be simulated
end
if isempty(nComp)
    nComp = 40;   % number of RGCs to be included in a computational run
end
if isempty(nRuns)
    nRuns = 1;     % number of computational runs
end
pRet.Ce=eye(d(1)*d(2)); % covariance matrix for Gaussian-induced noise correlations
pRet.magnitude = 1;      % parameters governing the nonlinearity mapping
pRet.gain      = 1;      % from linear filter responses to RGC output
pRet.offset    = -2.944;  % spiking probability
mode.tiling = 'random';           % arrangement of RGC centres 
mode.RF = 'center-surround DoG';  % layout of RGC receptive fields


idxC = zeros(nComp, nRuns);
for i = 1:nRuns
    idxC(:,i) = randsample(n,nComp)'; % index of nComp cells chosen for a run
end


% Settings concerning input images to be fed into the model retina
%--------------------------------------------------------------------------
if isempty(N)
    N = 10000; % number of image frames to be generated
end
if isempty(Nc)
    Nc = 500;  % chunks of images to be generated at a time (memory...)
end
if isempty(alpha)
    alpha = 2; % 0 'should' be white, 1 pink and 2 brown noise 
end
if isempty(xMag)
    xMag  = 1000; % magnitude of Gaussian noise to generate images from
end

% Settings concerning the maxEnt model fitting and evaluation
%--------------------------------------------------------------------------
modelFit = 'ising_count_l_0'; % maxEnt model to be fit to data
if isempty(nChains)
    nChains = 1; % number of individual MCMC chains to evaluate fitting results
end
if isempty(nSamplesEval)
    nSamplesEval = 20000; % number of Gibbs samples to extimate means E[f(X)]
end
if isempty(burnIn)
    burnIn       =  10000; % number of first Gibbs samples to be discarded
end
thinning     =     1; % distance in sequence between Gibbs samples to be
                      % stored (integers >1 thin out the MCMC chain)
                      
                      
fitoptions = struct;
fitoptions.optTol=1e-100; 
fitoptions.progTol=1e-100; 
%fitoptions.display='off';
fitoptions.MaxIter=3000;
fitoptions.MaxFunEvals=fitoptions.MaxIter;

MLhackData   = true;

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

data.d = d; data.n = n; data.nComp = nComp; data.nRuns = nRuns;
data.Nc = Nc; data.pRet = pRet; data.mode = mode; data.idxC = idxC;
data.alpha = alpha; data.xMag = xMag;
data.pars = pars; data.idxON = idxON; data.idxOFF = idxOFF;
data.Sigma = Sigma; data.SigmaDFs = SigmaDFs;


%%------------------------------------------------------------------------
% (2) Generate input images and do the retina simulation
%%-------------------------------------------------------------------------
disp('Generating RGC filters and spiking data')

[data.W, data.RGCcen] = genFilters(d,n,mode,pars); 
data.W=sparse(data.W);

data.spikes = zeros(n,N);
for i = 1:floor(N/Nc)
  disp([' - chunk ', num2str(i), ' out of ', num2str(floor(N/Nc))])
  x = xMag * spatialPattern([d(1),d(2),Nc], -alpha);
  tmp = retSim(x,data.W,pRet);
  data.spikes(:,(i-1)*Nc+1:i*Nc) = tmp.spikes;
end
clear tmp x;
disp('Computing correlations')
data.spkCorrs = full(corr(data.spikes'));
data.spkCov   = full(cov(data.spikes'));
data.RFoverlap = full(data.W*data.W');

tmp = pwd;
disp('Saving data')
switch tmp(1:3)
    case 'C:\'   
      save(['C:\Users\Loki\Desktop\Criticality\Results\run_', num2str(nComp), '_', num2str(ceil(now)),'.mat'],...
     'data', 'results');
    otherwise 
      save(['home\marcel\criticality\results\run_', num2str(nComp), '_', num2str(now),'.mat'],...
     'data', 'results');
end        

%% ------------------------------------------------------------------------
% (3) Fit the statistical model to data and evaluate the quality of fit
%%-------------------------------------------------------------------------
disp('Fitting maxEnt model')
datac = data.spikes(idxC,:);
disp('- augmenting data for ML')
if MLhackData % current method to avoid infinities in the ML fitting
 datac(:,end+1:end+nComp) = triu(ones(nComp,nComp));
 datac(:,end+1)       = zeros(nComp,1);
end
disp('- computing initialization for minFunc')
% for comparison: evaluate independent model, i.e. J = 0,L = 0
EX = full(mean(datac,2));
idxL = nComp*(nComp+1)/2 + (1:nComp+1);
lambdaInd = zeros(nComp*(nComp+3)/2+1,1); 
mfxInd = zeros(size(lambdaInd));
lambdaInd(1:nComp) = log(EX./(1-EX));
lambdaInd(lambdaInd==Inf) =   1000; % fairly hacky solution in case
lambdaInd(lambdaInd==-Inf) = -1000; % EX(k) = 0 resp. EX(k) = 1

fitoptions.lambda0 = -lambdaInd;
results.lambdaInd = lambdaInd;
mfxInd(1:nComp) = EX;
tmp = histc(sum(bsxfun(@lt, EX, rand([nComp,N])),1),0:nComp+1)/N;
if ~isempty(idxL)
  mfxInd(idxL) =tmp(1:nComp+1);
end
clear EX tmp idxL tmp

disp('- data size is ')
size(datac)

datac(datac>1) = 1;
disp('- starting MPF fitting')
[lambdaMPF,~,~,~,~] = fit_maxent_mpf(datac',fitoptions);
lambdaMPF = -lambdaMPF; % conventions...

results.fitoptions = fitoptions; 
results.lambda = lambdaMPF; 

tmp = pwd;
disp('Saving parameter fit')
switch tmp(1:3)
    case 'C:\'   
      save(['C:\Users\Loki\Desktop\Criticality\Results\run_', num2str(nComp), '_', num2str(ceil(now)),'.mat'],...
     'data', 'results');
    otherwise 
      save(['home\marcel\criticality\results\run_', num2str(nComp), '_', num2str(now),'.mat'],...
     'data', 'results');
end        
%%
disp('Generating data from fitted model for evaluation')
mfxTrain = zeros(nComp*(nComp+3)/2+1,1);
for i = 1:floor(N/Nc) 
  [fxTrain, ~] = setup_features_maxent(datac(:, (i-1)*Nc+1:i*Nc)', modelFit);
   mfxTrain = mfxTrain + mean(fxTrain, 1)'; 
end
mfxTrain = mfxTrain / floor(N/Nc);
mfxEval = zeros(nComp*(nComp+3)/2+1,1); 
tmpPwd = pwd;
for i=1:nChains
    disp(['- Chain #', num2str(i), ' of ', num2str(nChains)]);
    switch tmpPwd(1:3)
     case 'C:\'
        tmp = maxEnt_gibbs_pair_C(nSamplesEval, burnIn, lambdaMPF, logical(randi(2,[nComp,1])-1), 'laptop');
     otherwise
        tmp = maxEnt_gibbs_pair_C(nSamplesEval, burnIn, lambdaMPF, nComp, 'cluster');
    end
 storeTmp(:,i) = tmp;
 if mean((tmp-mfxTrain).^2)<mean((mfxEval-mfxTrain).^2) % if better than current chain results
   mfxEval = tmp;
 end 
end
clear fxTrain; 

results.mfxEval = mfxEval;
results.mfxTrain = mfxTrain;

disp('Saving sampled feature means')
switch tmpPwd(1:3)
    case 'C:\'   
      save(['C:\Users\Loki\Desktop\Criticality\Results\run_', num2str(nComp), '_', num2str(ceil(now)),'.mat'],...
     'data', 'results');
    otherwise 
      save(['/home/marcel/criticality/results/run_', num2str(nComp), '_', num2str(ceil(now)),'.mat'],...
     'data', 'results');
end        

disp('Plotting results')
figure; plot(mfxTrain, 'bo-'); hold on; plot(storeTmp, 'k.-'); plot(mfxEval, 'co-')

mainExperiment([], [], [], 1, 150000, 1000, 2, 1000, 3, 25000, 10000)