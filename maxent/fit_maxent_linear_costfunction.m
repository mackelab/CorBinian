function    [L,dL,ddL]=fit_maxent_linear_costfunction(lambda,x,means,weights,penalties)
%objective function (likelihood) for fitting a discrete maximum entropy (or
%log-linear) model of the form P(x)=1/Z exp(lambda'*x) to means "means"

%input: parameters lambda
%matrix of all possible states x of the system
%means: means that we want to match
%weights: set of fixed weights to be added
%penalties: weights for penalizing squared values of lambda

%output: (normalized) negative log-likelihood of data with means "means", 
%and its gradient and hessian

if nargin<=4
    penalties=zeros(size(lambda));
end

[logP,logZ,P]=logPMaxEnt(x,lambda,[],weights);

L=logZ-means*lambda-mean(weights)+ 0.5*(penalties.*lambda)'*lambda;


%keyboard
dL= x'*P-means';
dL= dL+ penalties.*lambda;

if nargout==3
ddL=x'*bsxfun(@times,x,P);
ddL=ddL+ diag(penalties);
end
