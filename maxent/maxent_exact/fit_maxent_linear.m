function [lambda,logZ, logP, fitmeans,output]=fit_maxent_linear(x, means, fitoptions,weights,penalties)
%Finds the parameters of a maximum entropy model of the form
%P(x)=1/Z exp(\sum_i lambda_i x_i), such that the means of x under this model
%match some supplied means.
%
%
%inputs:
%x: N by d matrix of all possible feature values
%means: d by 1 vector of feature-means that we want to match
%fitoptions: options of fit, as specified in the function minFunc by Mark
%Schmidt
% weights:. Optional, a vector of weights (of the same length) as x to be
% added to the log- probabilities, i.e. that the probabilities are then
%P(x)=1/Z exp(\sum_i lambda_i x_i+ weight)
%penalites: Optional, a vector the same size as means, which indicates the 
%weights with wich each lambda entry is penalized. The penalty term is
% Q= 0.5*sum_k lambda(k)^2 penalites(k). If we want to interpret this as a
% Bayesian Gaussian prior with variance sigma^2(k), then 
%penalties(k)= sigma^{-2}(k)/N_samples, where N_samples is the number of data-points.  
%
%outputs:
%
%lambda: vector of parameters of Max-Ent model
%logZ: log-partition function of model
%fitmeans: fitted means
%output: output structure returned by the function minimizer minFunc()
%
%uses minFunc by Mark Schmidt

[N,d]=size(x);

if nargin==2
    fitoptions=[];
end

if nargin<=3
    weights=0;
end

if nargin<=4
    penalties=ones(d,1)*1e-8;
elseif numel(penalties)==1
    penalties=penalties*ones(d,1);
end

%if starting point is provided, use it:
if isfield(fitoptions,'lambda0')
    lambda=fitoptions.lambda0(:);
else
    %otherwise use zeros:
    lambda=zeros(d,1);
end





%somewhat stupid option to restart algorithm
if ~isfield(fitoptions,'restarts')
    fitoptions.restarts=1;
end

funObj=@(lambda)(fit_maxent_linear_costfunction(lambda,x,means,weights,penalties));

for k=1:numel(fitoptions.restarts)
    try
        warning off
        [lambda,f,exitflag,output] = minFunc(funObj,lambda,fitoptions);
        warning on
    catch
        keyboard
        lambda=lambda*nan;
        f=nan;
        exitflag=nan;
        output=struct;
        output.iterations=nan;
    end
end

output.fs=f;
output.exitflag=exitflag;


%return log probabilities, log partition function, probabilities, and
%fitted means:
[logP,logZ,P,fitmeans]=logPMaxEnt(x,lambda,[],weights);

