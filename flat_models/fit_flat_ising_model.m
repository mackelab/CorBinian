function [h,J,count_distrib,model]=fit_flat_ising_model(mu,rho,N,maxN, fitoptions);
%function [h,J,count_distrib,model]=fit_flag_ising_model(mu,rho,N,maxN);
%
%fit flat ising model to population data with given mean and correlation
%functional form is P(k)=1/Z (n choose k) exp(hk+Jk^2); 
%
%
%
[mucount,varcount]=meancorr_2_meanvar_count(mu,rho,N);
%
if nargin<=3
    maxN=N;
end
if nargin<=4
    fitoptions=[];
end


states=[0:maxN];
x=[states',states'.^2];
means=[mucount,varcount+mucount^2];

lognchoosek=gammaln(N+1)-gammaln((0:maxN)'+1)-gammaln(N-(0:maxN)'+1);

weights=lognchoosek;

[lambda,logZ, logP, fitmeans,output]=fit_maxent_linear(x, means, fitoptions,weights);

h=lambda(1);
J=lambda(2);
count_distrib=exp(logP)';

if nargout==4
    model.h=h;
    model.J=J;
    model.lambda=lambda;
    model.features=x;
    model.count_distrib=count_distrib;
    model.weights=weights;
    model.fit=output;
    model.means=means;
    model.fitmeans=fitmeans;
    model.logZ=logZ;
    model.logP=logP;
end
    
