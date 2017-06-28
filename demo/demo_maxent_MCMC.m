%Demo Script:  Algorithms for fitting Ising models to large neural
%populations with iterative scaling and MCMC-sampling. Code by Jakob Macke,
%Jakob.Macke@gmail.com, using some methods originally by T.Broderick et al.
%
%This code is neither extensively tested, nor commented. Please do not
%distribute.
%
%Our parametrization of the Ising-model probabilities is that
%P(x)=1/Z*exp(h'*x+0.5*x'*J*x), where each x is a multivariate binary
%vector with states 0 and 1, and J is a strictly upper-triangular matrix
%
%set up paths etc:
clear all
clc


%number of neurons, must be <=12 for this demo-script to work
N=10;
%(Most functions should also work on much larger neural populations, only
%the function CompZ does not)

%%generate some random couplings and biases
randn('seed',5);
rand('seed',5);
p.htrue= randn(N,1)-2;
p.Jtrue = 3*(rand(N)-.5);
p.Jtrue = triu(p.Jtrue,1);
% calculate partition function of true model, mean mutrue, covariance
% Sigmatrue, and the probability of each state under the true model
[Ztrue, PIsingtrue, states, mutrue, Sigmatrue] = CompZ(p.htrue,p.Jtrue);


%%now, set parameters for fitting options

%size of each MC-sample
p.MC=20000;
%maximal number of outer-loop iterations
p.outerls=1000;
p.innerls=10;
%maximal run-time of the algorithm: 5 minutes
p.maxtime=5*60;

%desired accuracy of the means, correlations, covariances
%(metrics:
%1) mean absolute difference in means (of MC sample vs input)
%2) mean absolute difference of off-diagonal 2-point correlations
%3) mean absolute difference of off-diagonal pairwise covariances
p.exit_conditions=ApproxFinishLines(p.MC,mutrue,Sigmatrue);

%the algorithm exists if all three criteria are met, or if the maximal
%runtime or the maximal number of outer loops are finished. The algorithm
%only tries exiting at the end of an outer loop.

%delete all temporary files on completion:
p.remove_files=1;
%save parameters after each outer loop, not only at the very end
p.collect_all_params=1;

%do the actual training: function requires mean and covariance of the data,
%and the fitting options in the struct p.
[results p]= TrainIsing(mutrue,Sigmatrue,p);

%calculate partition function etc on the resulting data:
[Zest, PIsing, states, muest, Sigest] = CompZ(p.htrue,p.Jtrue);



%%
%display results and compare to true model
MCstates=logical(results.MC);

[PMC,PNoHit]=CountStates(MCstates,states);

PMC2=qIsing(states,results.h,results.J);
Zest=sum(PMC2);
PMC2=PMC2/Zest;
[PsMarg,PsInd]=CalcMarginals(states,PIsing);
%close all

%%
close all
subplot(2,5,1)
plot(p.htrue,'r')
hold on
plot(results.h,'b')
title(['True vs estimated bias terms'])
legend('True','Est')

subplot(2,5,6)
plot(mutrue,'r')
hold on
plot(muest,'b');
plot(results.mu,'g')

title(['True, estimated and sampled means'])
legend('True','estimated','sampled')

subplot(2,5,2)
imagesc(p.Jtrue,[-3,3])
title(['True coupling matrix'])

subplot(2,5,3)
imagesc(results.J,[-3,3])
title('Estimated coupling matrix')

subplot(2,5,7)
imagesc(results.cov_true-diag(diag(results.cov_true)),[-.05,.05])
title(['True covariance'])

subplot(2,5,8)
imagesc(results.cov-diag(diag(results.cov)),[-.05,.05])
title('Estimated covariance')


subplot(2,5,[4,5])
plot(log(PIsing),log(PsInd),'k.')
hold on
plot(log(PIsing),log(PMC2),'g.')
plot(log(sort(PIsing)),log(sort(PIsing)),'b-')
xlabel('True probabilities')
ylabel('Estimated probabilities')
legend('Independent model','Ising model','location','northwest')
title('True vs. estimated Ising Probabilities')


subplot(2,5,[9,10])
semilogy((1:length(results.track.metrics)),results.track.metrics,'.-')
hold on
semilogy(repmat(p.exit_conditions',1,length(results.track.metrics))')
xlim([1,length(results.track.metrics)])
xlabel('Number of outer loops')
ylabel('Performance metrics')
title('Convergence of the algorithm')



