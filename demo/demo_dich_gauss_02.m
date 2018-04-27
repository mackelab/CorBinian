% Demo script for software from the paper
% 'Generating spike-trains with specified correlations',
%  Macke et al., Neural Computation 2008
%
% Instructions:
%  - Change to the code directory and run demo_dich_gauss.m
%  - The script will take you through the functions .
%  - Read along in the demo.m file to follow what's happening.
%

fprintf('\nDemo program for GENERATING SPIKE TRAINS WITH SPECIFIED SPATIO-TEMPORAL CORRELATIONS, Macke et al.\n')

T = 1e5;  % number of time bins
N = 3;    % data dimensionality

rng(42) % fix seed

%% Step 1: Create some spatio-temporally correlated binary toy data 
% somewhat silly example to generate target spatio-temporal correlations

% create some temporal filters (to give temporal correlations), 
h = { [-0.2, -0.1, 0.,   0.,   0.,   0.2,  0.1]', ...
      [-0.5,  0.5, 0.66, 0.33, 0.22, 0.11, 0.05]', ...
      [-0.05, 0.1, 0.1, -0.1, -0.1, -0.1,  0.0]'};
for i = 1:N
    h{i} = flipud(h{i}) - mean(h{i});
    h{i} = h{i} / sqrt(sum(h{i}.^2));
end
tau = length(h{1});

% create some additive correlated noise (to give spatial correlations)
R = [0.3, 0.2, -0.1; 0.2, 0.9, 0. ; -0.1, 0., 0.4];
x = randn(T,N) * chol(R); 

% simulate data
y = zeros(T,N); y(1:tau, :) = x(1:tau, :);
for t = tau+1:T
    for j = 1:N
        y(t,j) = h{j}' * y(t-tau:t-1,j) + x(t,j);
    end
end

% binarize through a threshold (to give sensible means)
S = y  > 1.0;

clearvars -except S
pause(1)


%% Step 2: Sample from spatio-temporally correlated Dichotomized Gaussian

K = 10; % maximum temporal lag
[T, N] = size(S);

% construct Hankel
[mu, Sigma] = calc_hankel_from_data(S,K);

% obtain \gamma, \Lambda with core DichGauss code 
[St,g,L] = sampleDichGauss02(mu,Sigma,T,N,K);   

pause(1)

%% Step 3: verify correct means and (spatio-temporal) correlations
mu    = mean(S ,1); % estimate mean
muHat = mean(St,1); % 

fprintf('\nGenerating spatia-temporally correlated binary variables (section 2.3)\n\n')
fprintf('Target mean X1:  %.2f     Estimated mean X1:  %.3f\n',mu(1),muHat(1))
fprintf('Target mean X2:  %.2f     Estimated mean X2:  %.3f\n',mu(2),muHat(2))
fprintf('Target mean X3:  %.2f     Estimated mean X3:  %.3f\n',mu(3),muHat(3))

figure('pos',[10 10 900 600])
subplot(2,2,1)
imagesc(cov(S))
title('empirical covariance matrix')
colorbar
subplot(2,2,2)
imagesc(cov(St))
title('DG covariance matrix')
colorbar

subplot(2,2,3)
xcorrs_train = zeros(K-1, 2); xcorrs_test =  zeros(K-1, 2);
for i = 1:size(xcorrs_train,1)
    for j = 1:N
        xcorrs_test(i,j) = corr(St(1:end-2*K,j), St(i+1:end-2*K+i,j));
        xcorrs_train(i,j) = corr(S(1:end-2*K,j), S(i+1:end-2*K+i,j));
    end
end
clrs = {'b', 'k', 'r'};
for j = 1:N
    plot(xcorrs_train(:,j), 'o-', 'color', clrs{j}); 
    hold on; 
    plot(xcorrs_test(:,j), 'x--', 'color', clrs{j});
end
title('time-lagged correlations (one trace per variable)')
xlabel('time-lag \tau')
clear i j clrs
