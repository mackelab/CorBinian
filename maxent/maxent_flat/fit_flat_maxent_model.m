function [lambda,count_distrib,model]=fit_flat_maxent_model(count_distrib)
%given a count distribution, fit the 'collective ' model as I call it, defined in
%the Tkacik et al 'the simplest maximum entropy model for collective'
%paper. We write this model as 
%P(x)= 1/Z exp( lambda' F(x)) where the feature function F(x) is N
%dimensional, and the k-th elemnt is 1 i x has k ones, and 0 otherwise. 
%Following the sensible convention that also Tkacik uses, we shoose that P(K=0) has 'zero energy',
%i.e. P(silence)= 1/Z; 

count_distrib=count_distrib(:);
N=numel(count_distrib)-1;
Pzero=count_distrib(1);
logZ=-log(Pzero);

lognchoosek=gammaln(N+1)-gammaln((0:N)'+1)-gammaln(N-(0:N)'+1);

lambda=logZ-lognchoosek(2:end)+ log(count_distrib(2:end));

count_distrib=exp([0;lambda]'+lognchoosek'-logZ);



if nargout==3
    model=flat_model_calc_stats(count_distrib);
    model.lambda=lambda;
    model.count_distrib=count_distrib;
    model.logZ=logZ;
end
    