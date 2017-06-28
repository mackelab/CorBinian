n = 50; 
idxRepet = 1;
i = 1;  

   sig2_l2 = 400;
   sig2_sm = 10;

 a = 800; tau = 250;
 fitoptions = struct;
 fitoptions.burnIn  = 200;
 fitoptions.maxIter = min( 10*(n*(n+3)/2), 1000);
 fitoptions.nSamples = [0;round( a * 2.^((1:fitoptions.maxIter)' / tau ))];
 fitoptions.nSamples = floor(fitoptions.nSamples/100)*100;
 fitoptions.maxInnerIter = 1;    
 fitoptions.lambda0 = NaN;   % catch this later 
 fitoptions.regular = 'l1';  % not that l1 - regularization currently doesn't include the V(K)-Terms
 fitoptions.machine  = 'laptop';
 fitoptions.nRestart = 1;
 fitoptions.modelFit = 'ising_count_l_0';


 beta = 0.0001*ones( n*(n+3)/2 + 1,1); % very weak regularization
 eps = [0.01;0.05;0.01]; 
 hJV = ones(3,1); % default: all parameters, i.e. full model
 ifbwVK = true;
 ifSave = false; 
 fname  = [];
 
 load('C:\Users\Loki\Desktop\Criticality\code_critical_retina\Heat curves\EfxFFFData.mat')
 EfxNow = Efx{n/10}(:,idxRepet); % discard all other runs
  fitoptions.lambda0 = zeros(n*(n+3)/2+1,1); 
  fitoptions.lambda0(1:n) = log(EfxNow(1:n)./(1-EfxNow(1:n)));
  fitoptions.lambda0(fitoptions.lambda0==Inf) =  1000; % fairly hacky solu-
  fitoptions.lambda0(fitoptions.lambda0==-Inf)= -1000; % tion if Efx(i) = 0

i=1;
 [lambdaHat(:,i), fD{i}] = iterScaling(EfxNow, fitoptions, beta, eps, ...
                                        [], ...
                                        ifSave, [1;0;1], ifbwVK, ...
                                        sig2_l2, sig2_sm);
