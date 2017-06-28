
function testVarEnergyVariance(n, idxRepet, nSamplesMax, nSamplesBatch, burnIn)
% This serves to systematically test the length of MCMC samples needed to
% achieve stable estimates of the varience of energy.
% might want to discard the very first estimate Efy(:,1), Moe(:,1), as they
% might have been obtained without any burnin.

if nargin < 2
    idxRepet = 1;
end
if nargin < 3
    nSamplesMax = 1000000000;
end
if nargin < 4
    nSamplesBatch = 1000;
end
if nargin < 5
    burnIn = 0;
end

ticTime = now;
 load('/home/marcel/criticalityIterScaling/results/final/lambdaHatSim.mat') % got to transfer a couple of lambdas to the cluster,
  % gives 1x10 cell 'lambdaHatSim'               % ideally those fit to the full model (h,J,V) so far
if idxRepet <= size(lambdaHatSim{n},2)   
 lambda = lambdaHatSim{n}(:,idxRepet);
else
 warning('\lambda for desired idxRepet not found - using first one found')
 tmp = lambdaHatSim{n}(:,1);
 for j = 2:size(lambdaHatSim{n},2)
  if max(abs(lambdaHatSim{n}(:,j))) > 0
    tmp = lambdaHatSim{n}(:,j);
  end
 end
 lambda = tmp; clear tmp
end

n = 10*n;  % move from file index to population size 

N = floor(nSamplesMax/nSamplesBatch);    %number of update batches                                        
Efy = zeros(n*(n+3)/2+1,N); 
x0  = zeros(n,N+1); 
MoE = zeros(2,N); % moments of energy 

EX = exp(lambda(1:n))./(1+exp(lambda(1:n))); % a maxEnt model with only h,
x0(:,1) = double(rand(n,1)<EX);                   % i.e. no parameters J, L

for i = 1:N
   tocTime = now;
   if (tocTime - ticTime) > 1/3
      break
   end
   disp([num2str(i),'/',num2str(N)])
   [Efy(:,i), MoE(:,i), x0(:,i+1)] = maxEnt_gibbs_pair_C(nSamplesBatch, burnIn, lambda, x0(:,i), 'cluster');
   if mod(i,20)==0
    save(['/home/marcel/criticalityIterScaling/results/energy/VarEn', num2str(n), ...
          'idxRepet', num2str(idxRepet), 'nSamplesMax', num2str(nSamplesMax), '.mat'], ...
          'Efy', 'MoE', 'x0', 'lambda', 'nSamplesBatch', 'burnIn')
   end
end
   x0 = x0(:,1:i); Efy = Efy(:,1:i); MoE = MoE(:,1:i);
    save(['/home/marcel/criticalityIterScaling/results/energy/VarEn', num2str(n), ...
          'idxRepet', num2str(idxRepet), 'nSamplesMax', num2str(nSamplesMax), '.mat'], ...
          'Efy', 'MoE', 'x0', 'lambda', 'nSamplesBatch', 'burnIn')
end
