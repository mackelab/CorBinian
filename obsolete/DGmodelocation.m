function [rmode, rmax, rmin]=DGmodelocation(mus, rhos, usegammalambda)
%mode of DG as function of firing rate and correlation (using latent gamma
%and lambda if third argument is supplied

if numel(mus)==1
    mus=repmat(mus, size(rhos));
elseif numel(rhos)==1
    rhos=repmat(rhos, size(mus));
end
if nargin==2 || usegammalambda==0
    for k=1:numel(mus)
        [gammas(k), lambdas(k)]=RhoToLambda(mus(k), rhos(k));
    end
else
    gammas=mus;
    lambdas=rhos;
end

gammas=gammas(:);
lambdas=lambdas(:);
mus=mus(:);
rhos=rhos(:);

%find mode
rmode=normcdf(gammas.*sqrt(1-lambdas)./(1-2*lambdas));

%identify min and max

rmax=rmode;
rmax(mus<=.5 & lambdas>=.5)=0;
rmax(mus>.5 & lambdas<=.5)=1;

rmin=rmode;
rmin(mus<=.5 & lambdas<=.5)=1;
rmin(mus>.5 & lambdas<=.5)=0;

%special case: symmetry
if max(abs(mus-.5)<1e-10)
    rmax=.5*double([lambdas<.5, lambdas<.5])+[0*(lambdas>.5), lambdas>.5];
    rmin=.5*double([lambdas>.5, lambdas>.5])+[0*(lambdas<.5), lambdas<.5];
end


