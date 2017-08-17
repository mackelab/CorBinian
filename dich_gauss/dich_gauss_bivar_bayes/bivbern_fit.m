function [marginals, joints] = bivbern_fit(bins, priorTheta, priorLambda, LUTp11)
%BIVBERN_FIT fit a dichotomized bivariate Gaussian process to a bivariate
%Bernoulli process.
%   bins        : a Tx2 matrix containing one spike train in each column
%   priorTheta  : prior distribution of the parameters theta
%   priorLambda : prior distribution of the parameter lambda
%   LUTp11      : lookup table containing the probability of observing two
%   coincident spikes, for every combination of the tree parameters theta1,
%   theta2 and lambda.

theta = priorTheta.sup_x; 
lambda = priorLambda.sup_x;
pdfTheta = priorTheta.pdf_x;
pdfLambda = priorLambda.pdf_x;

% 1) compute the likelihood
likf = bivbern_lik(bins, lambda, theta, LUTp11);

% 2) Create a 3D matrix for the prior joint distribution:
% dim 1) firing probability neuron 1
% dim 2) firing probability neuron 2
% dim 3) covariance of the latent gaussian process
priorJointPdf = bsxfun(@times,...
    pdfTheta * pdfTheta', reshape(pdfLambda, 1, 1, length(lambda)));

% 3) find the posterior multiplying the likelihood and the prior
postJointPdf = likf .* priorJointPdf;
% normalize the posterior
% 3.a) integrate out theta 1
tmp = trapz(theta, postJointPdf, 1);
tmp = squeeze(tmp);
% 3.b) integrate out theta 2
tmp = trapz(theta, tmp, 1);
% 3.c) integrate out thata 3
normf = trapz(lambda, tmp);
% divide by the integral
postJointPdf = postJointPdf / normf;

%% find the marginals
% probFiring = 1 - normcdf(theta);
% [~, IX] = sort(probFiring);
% probFiring = probFiring(IX);

marginals = cell(3, 1);
% posterior marginal of theta1
pdfTheta1 = squeeze(trapz(...
    lambda, trapz(theta, postJointPdf, 2), 3));
pdfTheta1(pdfTheta1 < 0) = 0;
pdfTheta1 = pdfTheta1 / trapz(theta, pdfTheta1);
marginals{1} = ProbDist(theta, pdfTheta1);
% posterior marginal of theta2
pdfTheta2 = squeeze(trapz(...
    lambda, trapz(theta, postJointPdf, 1), 3));
pdfTheta2(pdfTheta2 < 0) = 0;
pdfTheta2 = pdfTheta2 / trapz(theta, pdfTheta2);
marginals{2} = ProbDist(theta, pdfTheta2);
% posterior marginal of lambda
pdfLambda = squeeze(trapz(...
    theta, trapz(theta, postJointPdf, 1), 2));
pdfLambda(pdfLambda < 0) = 0;
pdfLambda = pdfLambda / trapz(lambda, pdfLambda);
marginals{3} = ProbDist(lambda, pdfLambda);

%% find the joint distributions
joints = cell(3,1);
% prob(theta1, lambda)
pdf = squeeze(trapz(theta, postJointPdf, 2));
joints{1} = ProbDist2D(theta, lambda, pdf);
% prob(theta2, lambda)
pdf = squeeze(trapz(theta, postJointPdf, 1));
joints{2} = ProbDist2D(theta, lambda, pdf);
% prob(theta1, theta2)
pdf = trapz(lambda, postJointPdf, 3);
joints{3} = ProbDist2D(theta, theta, pdf);

end
