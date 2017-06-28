% this code implements minimum probability flow learning for a fully connected Ising model.  see http://arxiv.org/abs/0906.4779

% Author: Jascha Sohl-Dickstein (2009)
% Web: http://redwood.berkeley.edu/wiki/Jascha_Sohl-Dickstein
% This software is made available under the Creative Commons
% Attribution-Noncommercial License.
% (http://creativecommons.org/licenses/by-nc/3.0/)

description = 'd=50, 10000 samples'

% initialize
d = 5; % number of units
nsamples = 1000; % number of training samples
maxlinesearch = 10000; % this number is excessive just to be safe!!!!!! learning works fine if this is just a few hundred
independent_steps = 10*d; % the number of gibbs sampling steps to take between samples

run_checkgrad = 0;

minf_options = [];
%options.display = 'none';
minf_options.maxFunEvals = maxlinesearch;
minf_options.maxIter = maxlinesearch;
if run_checkgrad
   d = 5;
   nsamples = 2;
   minf_options.DerivativeCheck = 'on';
else
   % make the weight matrix repeatable
   rand('twister',355672);
   randn('state',355672);
end

% choose a random coupling matrix to generate the test data
J = randn( d, d ) / sqrt(d) * 3.;
J = J + J';
J = J/2;
J = J - diag(diag(J)); % set the diagonal so all the units are 0 bias
J = J - diag(sum(J));
fprintf( 'Generating %d training samples\n', nsamples );
burnin = 100*d;
% and generate the test data ...
t_samp = tic();
Xall = sample_ising( J, nsamples, burnin, independent_steps );
t_samp = toc(t_samp);
fprintf( 'training sample generation in %f seconds \n', t_samp );

% randomly initialize the parameter matrix we're going to try to learn
% note that the bias units lie on the diagonal of J
Jnew = randn( d, d ) / sqrt(d) / 100;
Jnew = Jnew + Jnew';
Jnew = Jnew/2;

L = randn( d+1, 1) / d;

Lnew = L;


% perform parameter estimation
fprintf( '\nRunning minFunc for up to %d learning steps...\n', maxlinesearch );
t_min = tic();

%%%%%%%%%%% choose one of these two Ising model objective functions %%%%%%%%%%%%%
% K_dK_ising is slightly faster, and includes connectivity only to states
% which differ by a single bit flip.
%Jnew = minFunc( @K_dK_ising, Jnew(:), minf_options, Xall );
% K_dK_ising_allbitflipextension corresponds to a slightly modified choice
% of connectivity for MPF. This modified connectivity includes an
% additional connection between each state and the state which is reached
% by flipping all bits.  This connectivity pattern performs better in cases
% (like neural spike trains) where activity is extremely sparse.

lambda = [Jnew(:);Lnew(:)];

lambdanew = minFunc( @K_dK_ising_PK, lambda, minf_options, Xall );
JnewPK = reshape(lambdanew(1:d^2), size(J));
LnewPK = lambdanew(end-d:end);

Jnew = minFunc( @K_dK_ising, Jnew(:), minf_options, Xall );

Jnew = reshape(Jnew, size(J));
t_min = toc(t_min);
fprintf( 'parameter estimation in %f seconds \n', t_min );

fprintf( '\nGenerating samples using learned parameters for comparison...\n' );
Xnew = sample_ising( Jnew, nsamples, burnin, independent_steps );
fprintf( 'sample generation with learned parameters in %f seconds \n', t_samp );

% generate correlation matrices for the original and recovered coupling matrices
mns = mean( Xall, 2 );
Xt = Xall - mns(:, ones(1,nsamples));
sds = sqrt(mean( Xt.^2, 2 ));
Xt = Xt./sds(:, ones(1,nsamples));
Corr = Xt*Xt'/nsamples;
mns = mean( Xnew, 2 );
Xt = Xnew - mns(:, ones(1,nsamples));
sds = sqrt(mean( Xt.^2, 2 ));
Xt = Xt./sds(:, ones(1,nsamples));
Corrnew = Xt*Xt'/nsamples;

Jdiff = J - Jnew;
Corrdiff = Corr - Corrnew;
jmn = min( [J(:); Jnew(:); Jdiff(:)] );
jmx = max( [J(:); Jnew(:); Jdiff(:)] );
cmn = min( [Corr(:); Corrnew(:); Corrdiff(:)] );
cmx = max( [Corr(:); Corrnew(:); Corrdiff(:)] );

% show the original, recovered and differences in coupling matrices
figure();
subplot(2,4,1);
imagesc( J, [jmn, jmx] );
axis image;
colorbar;
title( '{J}_{ }' );
subplot(2,4,2);
imagesc( Jnew, [jmn, jmx] );
axis image;
colorbar;
title( '{J}_{new}' );
subplot(2,4,3);
imagesc( Jdiff, [jmn, jmx] );
axis image;
colorbar;
title( '{J}_{ } - {J}_{new}' );
subplot(2,4,4);
imagesc( JnewPK, [jmn, jmx] );
axis image;
colorbar;
title( '{J}_{new, V(K)}' );

% show the original, recovered and differences in correlation matrices
subplot(2,4,5);
imagesc( Corr, [cmn,cmx] );
axis image;
colorbar;
title( '{C}_{ }' );    
subplot(2,4,6);
imagesc( Corrnew, [cmn,cmx] );
axis image;
colorbar;
title( '{C}_{new}' );    
subplot(2,4,7);
imagesc( Corrdiff, [cmn,cmx] );
axis image;
colorbar;
title( '{C}_{ } - {C}_{new}' );    
subplot(2,4,8);
imagesc( Corrnew, [cmn,cmx] );
axis image;
colorbar;
title( '{C}_{new, V(K)}' );    

figure();
plot( Corr(:), Corrnew(:), '.' );
axis([cmn,cmx,cmn,cmx]);
axis image;
xlabel( '{C}_{ }' );
ylabel( '{C}_{new}' );
title( 'scatter plot of correlations for original and recovered parameters' );
