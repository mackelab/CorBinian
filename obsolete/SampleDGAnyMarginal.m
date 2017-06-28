function [samples,hists] = SampleDGAnyMarginal(gammas,Lambda,supports,Nsamples)

% [samples,hists]=SampleDGAnyMarginal(gammas,Lambda,supports,Nsamples)
%   Generate samples for a Multivariate Discretized Gaussian with parameters
%   "gammas" and "Lambda" and "supports". The number of samples generated is "Nsamples"
%
%   input and output arguments are as described in "DGAnyMarginal"
%
% Usage: 
%
% Code from the paper: 'Generating spike-trains with specified
% correlations', Macke et al., submitted  for publication
%
% www.kyb.mpg.de/bethgegroup/code/efficientsampling

d=size(Lambda,1);

if isempty(supports)
    for k=1:d
        supports{k}=[0:numel(gammas{k})-1];
    end
end
        
% check whether Lambda is admissable
[cc, p] = chol(Lambda);
if p > 0 
  warning(['Covariance matrix of the latent Gaussian has at least one negative eigenvalue. ' ...
            'Applying Higham-Correction (see help higham).'])
  Lambda = higham(Lambda,1e-10,1e5);
  cc = chol(Lambda+eye(size(Lambda))*1e-6);
end
cc=chol(Lambda);

B=randn(Nsamples,d)*cc;

for k=1:d
    [hists{k},dd]=histc(B(:,k),[-inf;gammas{k};inf]);
    hists{k}=hists{k}/Nsamples;
    samples(:,k)=supports{k}(dd);
    hists{k}=hists{k}(1:max(1,end-2));
end
