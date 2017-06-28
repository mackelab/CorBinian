% Demo script for software from the paper
% 'Generating spike-trains with specified correlations',
%  Macke et al., Neural Computation 2008
%  code originally by P Berens, heavily modified and changes to conventions by JHM 2014
%
% Instructions:
%  - Change to the code directory and run demo_dich_gauss.m
%  - The script will take you through the functions .
%  - After each step, the script will pause. 
%  - To continue, hit any button.
%  - Read along in the demo.m file to follow what's happening.
%

f = 1;
fprintf('\nDemo program for GENERATING SPIKE TRAINS WITH SPECIFIED CORRELATIONS, Macke et al.\n')





%% Step 1: Generating correlated Poisson variables via the method 
% in section 3.3 for positive correlations
%
% Generate samples from a correlated Poisson distribution with
% specified mean and covariance matrix. The covariance between the two
% distributions is positive and quite strong. The construction follows now
% section 3.3, so we truncate a Gaussian to obtain a distribution with
% Poisson marginals and a given covariance.

fprintf('\nStep 1: Generating positively correlated Poisson variables (section 3.3)\n\n')
fprintf('Computing...\n')

mu = [7,9]';      % set the mean
C = [7 3;3 9];    % set an admissable covariance matrix

[S,gamma,Lambda,joints2D] = DGPoisson(mu,C,1e5);   % generate sample via discretized Gaussian

% We obtain a set of samples S, the parameters of the hidden Gaussian
% variable (gamma and Lambda) and a structure containing the marginal
% distributions as well as the 2D joint distribution

% Next, we verify that the samples have desired mean and covariance
% structure. We see that this is indeed the case.

muHat = mean(S,2);
CHat = cov(S');


fprintf('Target mean X1:  %.2f     Estimated mean X1:  %.3f\n',mu(1),muHat(1))
fprintf('Target mean X2:  %.2f     Estimated mean X2:  %.3f\n',mu(2),muHat(2))
fprintf('Target cov X1X2: %.2f     Estimated cov X1X2: %.3f\n',C(1,2),CHat(1,2))

%disp('To plot the marginals and joint histogram, hit any key...')
pause(1)


figure(f)

% Now, we look at the marginal distribution we obtain from our sampling
% procedure. For comparison, we also plot the true 1D Poisson distributions
% as specified by the mean of the marginal (dotted). We see that they are
% very close to one another. 

subplot(2,2,1)
h1 = joints2D{1,1};   % marginal of X_1
h2 = joints2D{2,2};   % marginal of X_2
plot(0:length(h1)-1,h1,'-r','markersize',5), hold on
plot(0:length(h1),poisspdf(0:length(h1),mu(1)),'.r','markersize',20)
plot(0:length(h2)-1,h2,'-g','markersize',5)
plot(0:length(h1),poisspdf(0:length(h1),mu(2)),'.g','markersize',20)
t = legend('X_1','X_1 Poisson','X_2','X_2 Poisson');
set(t,'box','off')
axis square
xlabel('X'), ylabel('P(X)')
title(sprintf('Step 4: Correlated Poisson, Positive Correlation\nMarginal Distributions with Poisson\n distribution with correct mean'))

% Finally, we convince ourselves that the two variables are indeed
% positively correlated by looking at the 2D joint histogram. We can see
% that although the marginals are perfect Poissonian, most samples are
% concentrated along the main diagonal, indicating the positive
% correlation.

subplot(2,2,3)
hh = joints2D{2};   % joint histogram
imagesc(hh)
axis square
title('')
xlabel('X_1'), ylabel('X_2')

%disp('To proceed to negatively correlated Poisson variables, hit any key...')
%paus

%% Step 2: Generating correlated Poisson variables via the method 
% in section 3.3 for negative correlations
%
% Again, we generate samples from a correlated Poisson distribution with
% specified mean and covariance matrix. This time the covariance is
% negative though.

fprintf('\nStep 5: Generating negatively correlated Poisson variables (section 3.3)\n\n')
fprintf('Computing...\n')

mu = [7,9]';      % set the mean
C = [7 -3;-3 9];    % set an admissable covariance matrix

[S,gamma,Lambda,joints2D] = DGPoisson(mu,C,1e5);   % generate sample via discretized Gaussian

% Next, we verify that the samples have desired mean and covariance
% structure. We see that this is indeed the case. 

muHat = mean(S,2);
CHat = cov(S');

fprintf('Target mean X1:  %.2f     Estimated mean X1:  %.3f\n',mu(1),muHat(1))
fprintf('Target mean X2:  %.2f     Estimated mean X2:  %.3f\n',mu(2),muHat(2))
fprintf('Target cov X1X2: %.2f     Estimated cov X1X2: %.3f\n',C(1,2),CHat(1,2))

%disp('To plot the marginals and joint histogram, hit any key...')
pause(1)

figure(f)

% Again we compare the marginal distributions to the Poisson distribution
% with the same mean and we find them matching very well.

subplot(2,2,2)
h1 = joints2D{1,1};   % marginal of X_1
h2 = joints2D{2,2};   % marginal of X_2
plot(0:length(h1)-1,h1,'-r','markersize',5), hold on
plot(0:length(h1),poisspdf(0:length(h1),mu(1)),'.r','markersize',20)
plot(0:length(h2)-1,h2,'-g','markersize',5)
plot(0:length(h1),poisspdf(0:length(h1),mu(2)),'.g','markersize',20)
t = legend('X_1','X_1 Poisson','X_2','X_2 Poisson');
set(t,'box','off')
axis square
xlabel('X'), ylabel('P(X)')
title(sprintf('Step 5: Correlated Poisson, Negative Correlation\nMarginal Distributions with Poisson\n distribution with correct mean'))

% Nevertheless, the samples have a different structure as we can see from
% the 2D joint histogram. This time, the two variables are negatively
% correlated and when X_1 tends to take large values, X_2 tends to take low
% ones. 

subplot(2,2,4)
hh = joints2D{2};   % joint histogram
imagesc(hh)
axis square
title('')
xlabel('X_1'), ylabel('X_2')

f=f+1
%% Step 3: Generating over-dispersed variables via the method 
% in section 3.3
%
% This time, we specify the marginal distributions directly

fprintf('\nStep 5: Generating variables with arbitary marginals (section 3.3)\n\n')
fprintf('Computing...\n')

%specified support and pmfs of our RVs:
supports{1}=[0:20];
pmfs{1}=1./sqrt(10+supports{1});
pmfs{1}=pmfs{1}/sum(pmfs{1});
supports{2}=[0:10];
pmfs{2}=exp(-(supports{2}-5).^2/10);
pmfs{2}=pmfs{2}/sum(pmfs{2});
[meano(1),varo(1)]=calc_mean_var(pmfs{1},supports{1});
[meano(2),varo(2)]=calc_mean_var(pmfs{2},supports{2});
Sigma=[varo(1),3;3,varo(2)];


[S,gammas,Lambda,joints2D,hists] = DGAnyMarginal(pmfs,Sigma,supports,1000)


% Next, we verify that the samples have desired mean and covariance
% structure. We see that this is indeed the case. 

muHat = mean(S,1);
CHat = cov(S);

fprintf('Target mean X1:  %.2f     Estimated mean X1:  %.3f\n',mu(1),muHat(1))
fprintf('Target mean X2:  %.2f     Estimated mean X2:  %.3f\n',mu(2),muHat(2))
fprintf('Target cov X1X2: %.2f     Estimated cov X1X2: %.3f\n',C(1,2),CHat(1,2))

%disp('To plot the marginals and joint histogram, hit any key...')
pause(1)

figure(f)

% Again we compare the marginal distributions to the specified distribution
% and  we find them matching very well.

subplot(2,2,1)
h1 = joints2D{1,1};   % marginal of X_1
h2 = joints2D{2,2};   % marginal of X_2
plot(supports{1},h1,'-r','markersize',5), hold on
plot(supports{1},pmfs{1},'.r','markersize',20)
plot(supports{2},h2,'-g','markersize',5)
plot(supports{2},pmfs{2},'.g','markersize',20)
t = legend('X_1','X_1 Poisson','X_2','X_2 Poisson');
set(t,'box','off')
axis square
xlabel('X'), ylabel('P(X)')
title(sprintf('Step 5: Correlated Poisson, Negative Correlation\nMarginal Distributions with Poisson\n distribution with correct mean'))


subplot(2,2,3)
hh = joints2D{2};   % joint histogram
imagesc(hh)
axis square
title('')
xlabel('X_1'), ylabel('X_2')



