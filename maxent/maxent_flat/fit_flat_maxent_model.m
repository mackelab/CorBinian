function [lambda,count_distrib,model]=fit_flat_maxent_model(count_distrib)
%function [lambda,count_distrib,model]=fit_flat_maxent_model(count_distrib);
%
% fit flat maxent model to population data with given mean and correlation.
% functional form is P(x)= 1/Z exp( lambda' F(x)) where the feature function 
% F(x) is N-dimensional, and the k-th elemnt is 1 if x has k ones, and 0 otherwise. 
%
% see Tkacik et al 'the simplest maximum entropy model for collective'
% we follow the convention that P(K=0) has 'zero energy',
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
    
