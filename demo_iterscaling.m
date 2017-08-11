function testPackageCluster(fname,d,nSamplesTrain,nSamplesEval,burnIn, ...
                               nSamplesIterScaling, burnInIterScaling, ...
                               maxIter, maxInnerIter, beta, ...
                               lambda0, lambdaTrue, mfxTrain, ...
                               a, tau, ...
                               ifSave, ifVK, ifbwVK)
% Simulation setup
%--------------------------------------------------------------------------
folder = '/home/marcel/criticalityIterScaling/results/';
if nargin < 1 || isempty(fname)
 fnames = [folder, 'tmp.mat'];
else
 fnames = [folder, fname];
end
if nargin < 2 || isempty(d)
d=15; %simulate 15 dimensional problem
end
if nargin < 3 || isempty(nSamplesTrain)
nSamplesTrain = 10000; %[100,1000,10000]; %generate 1000 data-points;
end
if nargin < 4 || isempty(nSamplesEval)
nSamplesEval  = 10000; %[100,1000,10000,100000];
end
if nargin < 5 || isempty(burnIn)
burnIn        = 1000;  %[100,1000,10000,10000];
end
if nargin < 6 || isempty(nSamplesIterScaling)
 fitoptions.nSamples = 10000;    
else 
 fitoptions.nSamples = nSamplesIterScaling;
 clear nSamplesIterScaling;
end
if nargin < 7 || isempty(burnInIterScaling)
 fitoptions.burnIn  = 1000;
else
 fitoptions.burnIn = burnInIterScaling;
 clear burnInIterScaling;
end
if nargin < 8 || isempty(maxIter)
 fitoptions.maxIter = 10*d;
else
 fitoptions.maxIter = maxIter;
 clear maxIter
end
if nargin < 9 || isempty(maxInnerIter)
 fitoptions.maxInnerIter = 1;    
else
 fitoptions.maxInnerIter = maxInnerIter;
 clear maxInnerIter
end
if nargin < 10 || isempty(beta)
 beta = 0.00001*ones(d*(d+1)/2 + d+1,1);
end
if nargin < 11 || isempty(lambda0)
 useLambdaInd = true; 
else 
 useLambdaInd = false;
 fitoptions.lambda0 = lambda0;
 clear lambda0;
end
if nargin < 14 || isempty(a)
 a = fitoptions.nSamples;
end
if nargin < 15 || isempty(tau)
 tau  = Inf;
end
if nargin < 16 || isempty(ifSave)
    ifSave = true;
end
if nargin < 17 || isempty(ifVK)
    ifVK = true;
end
if nargin < 18 || isempty(ifbwVK)
    ifbwVK = true;
end
if nargin > 13 && ~isempty(a) && ~isempty(tau)
 fitoptions.nSamples = [0;round( a * 2.^((1:fitoptions.maxIter)' / tau ))];
 fitoptions.nSamples = floor(fitoptions.nSamples/100)*100;
end

% The following shall be fixed to always these value
% general sampling/data generation/evaluation related (IT overhead...)
modelTrue = 'ising_count_l_0'; % full model
modelFit  = 'ising_count_l_0'; % full model
trainMethod = 'iterativeScaling';
% iterative scaling-related
fitoptions.regular = 'l1';
fitoptions.machine  = 'cluster';
fitoptions.nRestart = 1;
fitoptions.modelFit = modelFit;

if nargin < 12 || isempty(lambdaTrue)
 h=0.25*randn(d,1)-3.5; %generate random bias terms;
 J= (3)*0.15*(randn(d)); J=triu(J,1); 
 lambda=hJ2lambda(h,J);
 switch modelTrue
     case 'ising_count_l_0'
         %L=randn(d+1,1);
          L=zeros(d+1,1);
          L(1:28) = [3.7858,1.3932,-0.1493,-1.0014,-1.2690,-1.2713, ...
                    -1.2733,-1.0631,-0.7464,-0.5362,-0.2728,-0.0624, ...
                     0.2010,0.4114,0.5686,0.9383,1.0423,1.1465,1.0380,...
                     1.0359,0.9275,0.7128,0.4983,-0.1944,-0.0903, ...
                    -0.0925,-1.4757,-0.5218];
          L(29:end) = -5; 
     case 'ising_count'
         L=randn(d,1);
     case 'ising'
         L = [];
 end
 L = L - L(1); % a rather recently discovered nasty extra degree of freedom
 lambdaTrue = [lambda;L];
end

if nargin < 13 || isempty(mfxTrain)
 newData = true;
else
 newData = false;
end

pars.d = d; 
pars.beta = beta; 
pars.nSamplesTrain = nSamplesTrain;
pars.nSamplesEval  = nSamplesEval;
pars.burnIn = burnIn;



compTimes = zeros(length(nSamplesTrain),2,6);

for r = 1:length(nSamplesTrain)

compTimes(r,1,:) = clock;


% generate training data
%--------------------------------------------------------------------------
  % Initialize training data-generating MCMC chain with a sample drawn from
  % a nested model (only h = lamdbdaTrue(1:d), i.e. no J, no L)
  EX = exp(lambdaTrue(1:d))./(1+exp(lambdaTrue(1:d))); % works ok for small
  x0 = double(rand(d,1)<EX);                           % entries of L,J

if newData % else we already have mfxTrain
  disp('Generating training data')
 [mfxTrain,~,~] = maxEnt_gibbs_pair_C(nSamplesTrain(r), burnIn(r), lambdaTrue, x0, fitoptions.machine);
 %xTrain = maxEnt_gibbs_pair(nSamplesTrain(r), burnIn(r), thinning(r), ...
 %                       lambdaTrue, x0, modelTrue, 'default', 'samples');
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


save(fnames, 'lambdaTrue', 'beta', ...
            'mfxTrain', ...
            'fitoptions', 'pars');

% train model
%--------------------------------------------------------------------------
disp('Fitting maxEnt model')
switch trainMethod
  case 'MPF'
    fitoptions = struct;
    fitoptions.optTol=1e-100; 
    fitoptions.progTol=1e-100; 
   %fitoptions.display='off';
    fitoptions.MaxIter=3000;
    fitoptions.lambda0 = -lambdaInd;
    [lambdaHat,logZ,logP,fitmeans,output] = fit_maxent_mpf(xTrain',fitoptions);
    lambdaHat = -lambdaHat;
  case 'iterativeScaling'
    
    %fitoptions.lambda0 = lambdaTrue;
    %fitoptions.lambda0(1:(d*(d+1)/2)) = 0;
    %fitoptions.lambda0(1:d) = lambdaInd(1:d);
    
    %fitoptions.lambda0 = lambdaTrue;
    %fitoptions.lambda0(end-d:end) = 0;
    
    [lambdaHat, fitDiagnostics] = iterScaling(mfxTrain, fitoptions, beta, fname, ifSave, ifVK, ifbwVK);
end
% validate model
%--------------------------------------------------------------------------
 disp('Generating data from model fit')
 [mfxEval,~,~] = maxEnt_gibbs_pair_C(nSamplesEval(r), burnIn(r), lambdaHat, x0, fitoptions.machine);
%--------------------------------------------------------------------------
% if d < 20
%  [features,description,x]=setup_features_maxent(d,modelTrue);
%  [~,~,Ptrue, ~]=logPMaxEnt(features,lambdaTrue);
%  EX = sum(bsxfun(@times, x', Ptrue'),2);
%  description(isnan(description)) = d+1;
%  x1 = x; 
%  x1(:,end+1) = 1; 
%  EXX = sum(bsxfun(@times, (x1(:,description(1,d+1:d*(d+1)/2))...
%                         .* x1(:,description(2,d+1:d*(d+1)/2)))',Ptrue'),2);
%  EK = zeros(length(L),1);
% switch modelTrue
%     case 'ising_count'
%      for k = 1:length(EK)
%       EK(k) = sum((sum(x,2)==(k)) .* Ptrue);
%      end
%     case 'ising_count_l_0'
%      for k = 1:length(EK)
%       EK(k) = sum((sum(x,2)==(k-1)) .* Ptrue);
%      end
% end
% mfxTrue = [EX(:);EXX(:);EK(:)];
% clear x1 EX EXX EK description features Ptrue 
% end % if d < 20  
                    
compTimes(r,2,:) = clock;

end % end for r = 1:length(nSamplesTrain)

if ~exist('mfxTrue', 'var')
    mfxTrue = [];
end

save(fnames, 'lambdaTrue', 'lambdaHat',          'lambdaInd', ...
             'mfxTrain',   'mfxEval', 'mfxTrue', 'beta', ...
             'fitoptions', 'fitDiagnostics', 'pars', 'compTimes');

end