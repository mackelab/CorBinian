
% Simulation setup
%--------------------------------------------------------------------------
folder = '.';
fname = 'tmp';
fnames = [folder, 'tmp.mat'];
d=15;                  % data dimensionality
nSamplesData  = 10000; % draw from ground-truth parameters
nSamplesEval  = 10000; % draw from paramter estimates for comparison
burnIn        = 1000;  
fitoptions.nSamples = 10000;    
fitoptions.burnIn  = 1000;
fitoptions.maxIter = 1000;
fitoptions.maxInnerIter = 1;    
useLambdaInd = true; 
a = fitoptions.nSamples;
tau  = Inf;
eps = []; % loads defaults
beta = 0.00001*ones(d*(d+1)/2 + d+1,1);
hJV = [true;true;true];
ifSave = false;
ifVK = true;
ifbwVK = true;
fitoptions.nSamples = [0;round( a * 2.^((1:fitoptions.maxIter)' / tau ))];
fitoptions.nSamples = floor(fitoptions.nSamples/100)*100;

modelTrue = 'ising_count_l_0'; % full K-pairwise model
modelFit  = 'ising_count_l_0'; % 
trainMethod = 'iterativeScaling';

% iterative scaling-related
fitoptions.regular = 'l1';
fitoptions.machine  = 'cluster';
fitoptions.nRestart = 1;
fitoptions.modelFit = modelFit;

h=0.25*randn(d,1)-3.5; %generate random bias terms;
J= (3)*0.15*(randn(d)); J=triu(J,1); 
lambda=hJ2lambda(h,J);
switch modelTrue
    case 'ising_count_l_0'
      V=zeros(d+1,1);
      V(1:28) = [3.7858,1.3932,-0.1493,-1.0014,-1.2690,-1.2713, ...
                -1.2733,-1.0631,-0.7464,-0.5362,-0.2728,-0.0624, ...
                 0.2010,0.4114,0.5686,0.9383,1.0423,1.1465,1.0380,...
                 1.0359,0.9275,0.7128,0.4983,-0.1944,-0.0903, ...
                -0.0925,-1.4757,-0.5218];
      V(29:end) = -5; 
      V = V(1:d+1);
    case 'ising_count'
     V=randn(d,1);
    case 'ising'
     V = [];
end
 V = V - V(1);
 lambdaTrue = [lambda;V];

newData = true;

pars.d = d; 
pars.beta = beta; 
pars.nSamplesData = nSamplesData;
pars.nSamplesEval  = nSamplesEval;
pars.burnIn = burnIn;

compTimes = zeros(length(nSamplesData),2,6);

for r = 1:length(nSamplesData)

compTimes(r,1,:) = clock;


% generate training data
%--------------------------------------------------------------------------
  % Initialize training data-generating MCMC chain with a sample drawn from
  % a nested model (only h = lamdbdaTrue(1:d), i.e. no J, no L)
  EX = exp(lambdaTrue(1:d))./(1+exp(lambdaTrue(1:d))); % works ok for small
  x0 = double(rand(d,1)<EX);                           % entries of L,J

  if newData % else we already have mfxTrain
  disp('Generating training data')
  [mfxTrain,~,~] = maxEnt_gibbs_pair_C(nSamplesData(r), burnIn(r), ...
                                       lambdaTrue, x0, fitoptions.machine);
  end

 EX = mfxTrain(1:d);
 
 lambdaInd = zeros(size(lambdaTrue)); 
 lambdaInd(1:d) = log(EX./(1-EX));
 lambdaInd(lambdaInd==Inf) =   1000; % fairly hacky solution in case
 lambdaInd(lambdaInd==-Inf) = -1000; % EX(k) = 0 resp. EX(k) = 1

 if useLambdaInd 
   fitoptions.lambda0 = lambdaInd;
 end
 clear EX tmp xTrain fxTrain 


% train model
%--------------------------------------------------------------------------
disp('Fitting maxEnt model')
switch trainMethod
  case 'MPF'
    fitoptions = struct;
    fitoptions.optTol=1e-100; 
    fitoptions.progTol=1e-100; 
    fitoptions.MaxIter=3000;
    fitoptions.lambda0 = -lambdaInd;
    [lambdaHat,logZ,logP,fitmeans,output] = fit_maxent_mpf(xTrain',fitoptions);
    lambdaHat = -lambdaHat;
    
  case 'iterativeScaling'
      
    [lambdaHat, fitDiagnostics] = iterScaling(mfxTrain, fitoptions, ...
                              beta, eps, fname, ifSave, hJV, ifVK, ifbwVK);
end


% validate model
%--------------------------------------------------------------------------
 disp('Generating data from model fit')
 [mfxEval,~,~] = maxEnt_gibbs_pair_C(nSamplesEval(r), burnIn(r), ...
                                     lambdaHat, x0, fitoptions.machine);

if d < 20
 [features,description,x]=setup_features_maxent(d,modelTrue);
 [~,~,Ptrue, ~]=logPMaxEnt(features,lambdaTrue);
 EX = sum(bsxfun(@times, x', Ptrue'),2);
 description(isnan(description)) = d+1;
 x1 = x; 
 x1(:,end+1) = 1; 
 EXX = sum(bsxfun(@times, (x1(:,description(1,d+1:d*(d+1)/2))...
                        .* x1(:,description(2,d+1:d*(d+1)/2)))',Ptrue'),2);
 EK = zeros(length(L),1);
switch modelTrue
    case 'ising_count'
     for k = 1:length(EK)
      EK(k) = sum((sum(x,2)==(k)) .* Ptrue);
     end
    case 'ising_count_l_0'
     for k = 1:length(EK)
      EK(k) = sum((sum(x,2)==(k-1)) .* Ptrue);
     end
end
mfxTrue = [EX(:);EXX(:);EK(:)];
clear x1 EX EXX EK description features Ptrue 
end % if d < 20  
                    
compTimes(r,2,:) = clock;

end % end for r = 1:length(nSamplesData)

if ~exist('mfxTrue', 'var')
    mfxTrue = [];
end
%%
figure; 
subplot(131)
plot(mfxTrue(1:d), mfxEval(1:d), '.');
hold on
plot(mfxTrain(1:d), mfxEval(1:d), 'k.');
title('first moments')
axis square
subplot(132)
plot(mfxTrue(d+1:end-d-1), mfxEval(d+1:end-d-1), '.')
hold on
plot(mfxTrain(d+1:end-d-1), mfxEval(d+1:end-d-1), 'k.')
title('second moments')
axis square
subplot(133)
plot(mfxTrain(end-d:end), 'ko-')
hold on
plot(mfxTrue(end-d:end), 'o-')
plot(mfxEval(end-d:end), 'ro-')
hold off
legend({'data', 'true', 'est.'})
title('population counts')
