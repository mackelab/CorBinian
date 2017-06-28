function [h,logZ,count_distrib,model]=flat_model_vary_temp_fix_mean(mean_rate,count_distrib,beta,fitoptions)
%given a count distribution, fit the 'collective ' model as I call it, defined in
%the Tkacik et al 'the simplest maximum entropy model for collective'
%paper. We write this model as
%P(x)= 1/Z exp(hk +alpha* lambda' F(x)),, where F(x) etc are defined as in
%fit_flat-maxent_model, and alpha and lambda are given. We waht to adjust Z
%such that the model has the right mean firing rates.

N=numel(count_distrib)-1;

if nargin<=3 fitoptions=[]; end


%we find the best fitting h by solving a maximum entropy problem with
%weights,i.e. fit P(K=k)=1/Z exp(h*k+ weight(k)) for fixed weights.
lognchoosek=(gammaln(N+1)-gammaln((0:N)+1)-gammaln(N-(0:N)+1))';
[lambda,count_distrib_check,model]=fit_flat_maxent_model(count_distrib);

lambda(lambda==-Inf) = -1000;

weights=beta*[0;lambda]+lognchoosek;

x=[0:N]';
means=mean_rate*N;
[h,logZ, logPK, fitmeans,output]=fit_maxent_linear(x, means, fitoptions,weights);

count_distrib=exp(logPK');
%keyboard

if nargout==4
    model=flat_model_calc_stats(count_distrib);
    model.h=h;
    model.logZ=logZ;
    model.fitmean=fitmeans;
end