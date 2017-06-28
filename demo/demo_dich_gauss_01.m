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


%% Step 1: Generating correlated binary variables (2D example)
mu = [.4,.3]';      % set the mean
v = mu.*(1-mu);     % calculate the variances
C = diag(v);        % build covariance matrix
C(1,2) = .1;
C(2,1) = .1;

[S,g,L] = sampleDichGauss01(mu,C,1e5);   % generate samples from the DG model

muHat = mean(S,1);  % estimate mean
CHat = cov(S);     % estimate covariance

fprintf('\nStep 1: Generating correlated binary variables (section 2.1)\n\n')
fprintf('Target mean X1:  %.2f     Estimated mean X1:  %.3f\n',mu(1),muHat(1))
fprintf('Target mean X2:  %.2f     Estimated mean X2:  %.3f\n',mu(2),muHat(2))
fprintf('Target cov X1X2: %.2f     Estimated cov X1X2: %.3f\n',C(1,2),CHat(1,2))

%disp('To proceed to compare histograms, hit any key...')
pause(1)

%% Step 2: How much do correlations change the distribution when means are
% assumed to be equal?

fprintf('\nStep 2: Effect of correlations\n\n')
fprintf('Computing...\n')

% First, we generate a mean vector in ten dimensions and a covariance
% matrix:

g=randn(10,1)-1;
G=randn(10);
Lambda=cov_2_corr(G*G');
[mu,C]=Lambda_2_Sigma(g,Lambda);



% We find the histogram of the independent distribution with this mean
% vector P(x) = PROD_i P(X_i=1) and the DG distribution with the same mean
% and the covariances as defined above.

states=all_states(10);
S = sampleDichGauss01(mu,C,1e5);
h1=CalcIndep(mu,states);
[h2]=CountStates(S,states);


% When we compare them, we see that some patterns occur much more often 
% in the correlated then in the independent distribution. This a clear
% indication of the strength of the correlations. The "stripes" coincide
% with patterns with a different number of active neurons (starting at 0 in
% the upper left corner).

figure(f)
loglog(h1,h2,'k.')
xlabel('P(X) in independent model')
ylabel('P(X) in DG model')
title('Step 2: What effect do correlations have?')
axis square

f = f+1;


