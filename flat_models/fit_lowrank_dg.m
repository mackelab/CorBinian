function [gamma,Lambda, Loadings, ResVar]=fit_lowrank_dg(mu,Sigma,d,acc);
%function [gamma,Lambda,Loadings,ResVar]=fit_lowrank_dg(mu,Sigma,d,acc);
%
%given mean firing rates mu and the covariance Sigma, fit a DG model to the
%data which is constrained to be low-rank of dimensionality d. acc is the
%accurary of the covariance calcuation (see findLatentGaussian01),
%tune_iter the number of 'fine-tuning' runs to do. In each fine-tuning run,
%the 'worst' covariance estimated is fine-tuned by hand.
%
% definitiion: x=1 iff z>0, where z ~ N(gamma,Lambda), Lambda of rank d and
% with 1 on the diagonal

%function is super-naive and inefficient at the moment--very easy to make
%it much faster by getting rid of some of the for-loops.

if nargin==3
    acc=1e-5;
end

    

%assert(d==1, 'Only implemented for 1-D models so far!!!')

%first, find full underlying covariance, i.e. assuming full rank
[gamma,Lambda] = findLatentGaussian01(mu,Sigma,acc);
%ensure that Lambda is positive definite:
Lambda_pd = higham(Lambda);


%then, apply factor analysis to Lambda to reduce it to one-d matrix
[Loadings, ResVar]=factoran(Lambda_pd,d,'xtype','covariance');

Lambda_fa=Loadings*Loadings'+diag(ResVar);



%renormalize to 1: 
vars=diag(Lambda_fa);
ResVar_norm=max(0,ResVar+(1-vars));
Lambda_fa_norm=Loadings*Loadings'+diag(ResVar_norm);


Lambda=Lambda_fa_norm;
ResVar=ResVar_norm;




%% fine-tuning not implemented yet!!!
%keyboard


%Lambda_tuned=Lambda_fa_norm;
%in the end, do a bit of fine-tuning:
%for k=1:tune_iter
    %calculate difference between original and reconstructed covariance:
%    [mu_test, Sigma_test]=Lambda_2_Sigma(gamma, Lambda_tuned);
%    delta=Sigma-Sigma_test;
    
%    delta=
%    if max(abs(delta))<acc
%       break 
%    end
%end


%Lambda=Lambda_tuned;



