function runMaxEntSimData(n, idxRepet, fitoptions, beta, eps, a,tau, hJV, fname, sig2_l2, sig2_sm)
% Fits a maximum entropy model with parameters for feature functions f(X)
% to data from the simulation of a population of RGCs. 
% Loads pre-computed data and applies subsets of the model 

n = 10*n;  % if only growing populations in steps of 10, this makes the 
           % relation between population size and cluster core much easier. 

% input:
%  n:        number of cells from the simulation to include in model fit
%  idxRepet: index of pre-computed set of data moments E[f(X)] to use
%            (can be of 50>=length(idxRepet)>=1, will fit models to several 
%             sets of data moments one after the other)
% fitoptions: (optional) structure specifying options for the model fits
% beta:       (optional) vector of regularization weights for the models
% a, tau:     (optional) values for computing the length of individual
%                        MCMC samples by nSamples(iter) = a*2^(iter/tau)
% hJV:        (optional) 3-by-1 vector of booleans giving whether the
%                        respective model parts h, J, V are to be included
% fname:      filename for stored results. Will automatically be appended
%             with number of cells n and number of data set from idxRepet

if nargin < 3 || isempty(fitoptions)
 fitoptions.nSamples = 3200; % may be overwritten below by 'a' and 'tau'  
 fitoptions.burnIn  = 200;
 fitoptions.maxIter = min( 100*(n*(n+3)/2), 10000);
 fitoptions.maxInnerIter = 1;    
 fitoptions.lambda0 = NaN;   % catch this later 
 fitoptions.regular = 'l1';  % not that l1 - regularization currently doesn't include the V(K)-Terms
 fitoptions.machine  = 'cluster';
 fitoptions.nRestart = 1;
 fitoptions.modelFit = 'ising_count_l_0';
 
if nargin < 4 || isempty(beta)
 beta = 0.0001*ones( n*(n+3)/2 + 1,1); % very weak regularization
end

if nargin < 5 || isempty(eps)
 eps = [0.01;0.05;0.01]; % tolerance for E[X] (h), cov(X) (J), P(K) (V(K)).  
end

if nargin > 6 && ~isempty(a) && ~isempty(tau)
 fitoptions.nSamples = [0;round( a * 2.^((1:fitoptions.maxIter)' / tau ))];
 fitoptions.nSamples = floor(fitoptions.nSamples/100)*100;
end
  
if nargin < 8 || isempty(hJV)
 hJV = ones(3,1); % default: all parameters, i.e. full model
end

folder = '/home/marcel/criticalityIterScaling/results/'; % output folder
if nargin < 9 || isempty(fname)    
 fname = date;           % note that fname will be further 
end                      % augmented by the iteration number
fname = [fname, 'n', num2str(n)]; 


if nargin < 10 || isempty(sig2_l2)
   sig2_l2 = 400;
end

if nargin < 11 || isempty(sig2_sm)
   sig2_sm = 100;
end


%load('/home/marcel/criticalityIterScaling/data/EfxCB2Data.mat')
%load('/home/marcel/criticalityIterScaling/data/EfxFFFData.mat')
load('/home/marcel/criticalityIterScaling/data/EfxSimData.mat')
% gives Efx (data), idxSubsample (neuron index table), par (for check-up)

ifSave = true; % always store results 
ifbwVK = true; % always do block-wise update 

lambdaHat = zeros(n*(n+3)/2+1,size(Efx,2));
fD        = cell(size(Efx,2),1); % might run out of memory for many repets

EfxNow = zeros(n*(n*3)/2+1, 1);
cd /home/marcel/criticalityIterScaling/results/
for i = idxRepet

 % find data moments E[f(X)] for currently chosen n and idxRepet
 for j = 1:length(idxSubsamples)
  if size(idxSubsamples{j},1)==n
      EfxNow = Efx{j}(:,i); % discard all other runs
      break
  end
 end

  fitoptions.lambda0 = zeros(n*(n+3)/2+1,1); 
  fitoptions.lambda0(1:n) = log(EfxNow(1:n)./(1-EfxNow(1:n)));
  fitoptions.lambda0(fitoptions.lambda0==Inf) =  1000; % fairly hacky solu-
  fitoptions.lambda0(fitoptions.lambda0==-Inf)= -1000; % tion if Efx(i) = 0
 
 disp(['fitting on data set', num2str(i), '/', num2str(size(Efx,2))])   
 mkdir([fname, 'idxRepet', num2str(i)])
 [~, ~] = iterScaling(EfxNow, fitoptions, beta, eps, ...
                                       [fname, 'idxRepet', num2str(i)], ...
                                       ifSave, hJV, ifbwVK, ...
                                       sig2_l2, sig2_sm);
 
end

end % end function
