clc
%clear all
%close all

% Set simulation up
%--------------------------------------------------------------------------
d=10; %simulate 10 dimensional problem
n= 1000; %generate 1000 data-points;
nSamples = 1000; burnIn = 1000; 
thinning = 1;
model = 'ising_count_l_0';
newLambda = true;
nRestarts = 25;

if newLambda
 h=randn(d,1)-1; %generate random bias terms;
 J=randn(d); J=triu(J,1)*10;%/sqrt(d); 
 %J = zeros(size(J));
  
 lambda=hJ2lambda(h,J);
 switch model
     case 'ising_count_l_0'
         L=randn(d+1,1)/sqrt(d);
     case 'ising_count'
         L=randn(d,1)/sqrt(d);
     case 'ising'
         L = [];
 end
 lambda = [lambda;L];
end
figure;  
subplot(2,16,1:16) % subplot for true and recovered moments E[f(X)]
title('A: True and recovered distribution moments E[f_i(X)]')
xlabel('i'); ylabel('$$E[f_i(X)]$$', 'Interpreter', 'Latex')  


% Pre-computations
%--------------------------------------------------------------------------
if d < 20
 [features,description,x]=setup_features_maxent(d,model);
 [logPtrue,logZtrue,Ptrue, means_true]=logPMaxEnt(features,lambda);
 EX = sum(bsxfun(@times, x', Ptrue'),2);
 description(isnan(description)) = d+1;
 x1 = x; x1(:,end+1) = 1; 
 EXX = sum(bsxfun(@times, (x1(:,description(1,d+1:d*(d+1)/2))...
                         .*x1(:,description(2,d+1:d*(d+1)/2)))',Ptrue'),2);
 EK = zeros(length(L),1);
 switch model
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
 clear x1 EX EXX EK

% Generate directly sampled data (via inverse CDF essentially) to compare
%--------------------------------------------------------------------------
 x_sampled=sample_discrete(x,Ptrue,max(n));
 fx_sampled = setup_features_maxent(x_sampled, model);
 fx_sampled = full(fx_sampled);
 %x0 = x_sampled(randi(n),:)';         % starting point for Gibbs sampler
 x0 = d;
 
% Plot pre-computation and discrete-sampling results
%--------------------------------------------------------------------------
 plot(mfxTrue,  'r--', 'linewidth', 1.5);
 hold on
 mfxTrain = mean(fx_sampled,1)'; 
 plot(mfxTrain, 'ko-', 'linewidth', 2)
 hold off

else % need a starting point for Gibbs sampler also in high dimensions...
  x0 = d; % tells maxEnt_gibbs to draw its own initial MCMC chain element

end % if d < 20, pre-computation and discrete-sampling part

% Start MCMC chains
%--------------------------------------------------------------------------
mfxEval = zeros(length(lambda), nRestarts);
for i = 1:nRestarts
disp(['Restart number ', num2str(i), ' of ', num2str(nRestarts)])
% Do the sampling (the part that takes a while...)
[xSampled] = maxEnt_gibbs(nSamples, burnIn, thinning, lambda, x0, model);
fxSampled = setup_features_maxent(xSampled', model);
fxSampled = full(fxSampled);

% Update figure:
hold on
mfxEval(:,i) = mean(fxSampled(end/2:end,:),1)'; % second half of data
plot( mfxEval(:,i), 'go-', 'linewidth', 1)
mfxEval(:,i) = mean(fxSampled,1)';
plot( mfxEval(:,i), 'bo-', 'linewidth', 2)
hold off
%clear fxSampled xSampled % huge memory consumption here
end
% Finish figure:
hold on
if d < 20
 plot(mfxTrain, 'ko-', 'linewidth', 2.5) % replot these one 
 plot(mfxTrue,  'r--', 'linewidth', 1.5) % for enhanced visibility  
 axis([0.5, size(fxSampled,2)+0.5, ...
           min([mfxTrain(:);mfxEval(:)]),max([mfxTrain(:);mfxEval(:)])])
 line([d,d]+            0.5,[min([mfxTrain(:);mfxEval(:)]), ... 
                 max([mfxTrain(:);mfxEval(:)])], 'color', 'k') 
 line([d,d]+(d*(d-1)/2)+0.5,[min([mfxTrain(:);mfxEval(:)]), ... 
                 max([mfxTrain(:);mfxEval(:)])], 'color', 'k') 

 legend('true', 'sampled discrete', 'second half of MCMC chain', ...
        'full MCMC chain', 'location', 'north') 
else
 axis([0.5, size(fxSampled,2)+0.5, min(mfxEval(:)),max(mfxEval(:))])
 line([d,d]+ 0.5,[min(mfxEval(:)),max(mfxEval(:))],'color','k')% delineate 
 line([d,d]+(d*(d-1)/2)+0.5,[min(mfxEval(:)), ...              % parameters
                             max(mfxEval(:))], 'color', 'k')   % for h,J,L
 legend('second half of MCMC chain','full MCMC chain', 'location', 'north') 
end
hold off
box off; set(gca, 'TickDir', 'out')
legend boxoff

% Sanity check: compare count distribution to one that is expected if one
% actually had fitted an independent model (J_ij = 0 for all i,j)
% means = mean( xSampled,2);
% tmp = zeros(nSamples,1);
% for t = 1:nSamples % for each samples ...
%  for k = 1:d % ... draw d times independently from E[X_k = 1] ...
%   tmp(t) = tmp(t) + (rand(1)<means(k)); % ... and add all up
%  end
% end
% plot(0:d, histc(tmp, 0:d)/nSamples, 'kx--');

% Add subplots for MSE comparison between different chains and 
% discretely sampled stuff. Always take 'true' moments as reference
if d<20

lbls = cell(nRestarts+1,1);
lbls{1} = 'Sampled';
for i = 1:nRestarts
 lbls{i+1} = ['#', num2str(i)];
end
subplot(2,16,21:24)
idx = 1:d;
mse = mean(bsxfun(@minus, mfxTrue(idx), ...
                         [mfxTrain(idx), mfxEval(idx,:)]).^2,1);
plot(1:nRestarts+1,mse, 'o--', 'linewidth', 2)
set(gca, 'XTick', 1:nRestarts+1);
set(gca, 'XTickLabel', lbls);
title('C: entries for h')
box off
axis([0.5, nRestarts+1.5, 0.9*min(mse), 1.1*max(mse)]) 
subplot(2,16,25:28)
idx = d+1:d*(d+1)/2;
mse = mean(bsxfun(@minus, mfxTrue(idx), ...
                         [mfxTrain(idx), mfxEval(idx,:)]).^2,1);
plot(1:nRestarts+1,mse, 'o--', 'linewidth', 2)
set(gca, 'XTick', 1:nRestarts+1);
set(gca, 'XTickLabel', lbls);
title('D: entries for J')
box off
axis([0.5, nRestarts+1.5, 0.9*min(mse), 1.1*max(mse)]) 
subplot(2,16,29:32)
idx = length(mfxTrue)-length(L)+1:length(mfxTrue);
mse = mean(bsxfun(@minus, mfxTrue(idx), ...
                         [mfxTrain(idx), mfxEval(idx,:)]).^2,1);
plot(1:nRestarts+1,mse, 'o--', 'linewidth', 2)
set(gca, 'XTick', 1:nRestarts+1);
set(gca, 'XTickLabel', lbls);
title('E: entries for L')
box off
axis([0.5, nRestarts+1.5, 0.9*min(mse), 1.1*max(mse)]) 
subplot(2,16,17:20)
idx = 1:length(mfxTrue);
mse = mean(bsxfun(@minus, mfxTrue(idx), ...
                         [mfxTrain(idx), mfxEval(idx,:)]).^2,1);
plot(1:nRestarts+1,mse, 'o--', 'linewidth', 2)
set(gca, 'XTick', 1:nRestarts+1);
set(gca, 'XTickLabel', lbls);
box off
axis([0.5, nRestarts+1.5, 0.9*min(mse), 1.1*max(mse)]) 
ylabel('$$\frac{1}{d} \Sigma_{k=1}^d (E[f_k(X)] - \hat{E}[f_k(X)])^2$$',...
       'Interpreter', 'Latex')
title('B: MSE for all entries of \lambda')

else 
subplot(2,16,17:24)
plot(mean(mfxEval,2), 'r--', 'linewidth', 2) 
hold on

end % if d < 20, plotting part 

% Finalize image by giving it a nice super-title
switch model
    case 'ising'
        model_name = 'Ising';
    case 'ising_count'
        model_name = 'Ising + V(K)';
    case 'ising_count_l_0'
        model_name = 'Ising + V(K)';
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', ['Gibbs sampling results for ', model_name, ', d = ', ...
                num2str(d), ', #Samples = ', num2str(nSamples), ...
                ', #burnIn = ', num2str(burnIn)], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
