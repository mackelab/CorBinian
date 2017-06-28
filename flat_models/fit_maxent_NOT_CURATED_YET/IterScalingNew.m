function [h,J,Zfinal,ps,f,hofx,beta,output]=IterScalingNew(mu,Sigma, n,acc,maxiter,typer,hJ0,Bx,EBx)

%hkeyboard
femp=[mu;(Sigma+mu^2)/2];

if nargin<=5
    typer='maxent';
end
if nargin<=4
    maxiter=500;
end
if nargin<=3 || isempty(acc)
    acc=10e-20;
end
if nargin>=8 && ~any(isnan(Bx))
    useBx=true;
else
    useBx=false;
end

if n>10000;
    h=nan;
    J=nan;
    Zfinal=nan;
    ps=nan;
    f=nan;
    hofx=nan;
    return
end
   


range=(0:n)';
featurevals=[range,range.^2/2]';
J=0;

switch typer
    case 'maxent' 
        p=1/(mu+1);
       h=0;
        hofx=0; % do not need to precompute combinatorial stuff

    case 'binom'
        p=mu/n;
        h=log(p/(1-p));
        J=Sigma/mu/(n-mu);
        J=1.5*J/n;
        hofx=(gammaln(n+1)-gammaln(range+1)-gammaln(n-range+1))';
end

if nargin>=7 && ~any(isnan(hJ0))
    h=hJ0(1);
    J=hJ0(2);
end
if numel(hJ0)==3
    beta=hJ0(3);
else
    beta=0;
end

lambda=[h;J];
if useBx
   lambda=[lambda; beta];
   featurevals=[featurevals;Bx(:)'];
   femp=[femp; EBx];
  % keyboard
end

%[lambda, f, nits] = Minimizer(lambda,@FitMaxEntObjFunNew,maxiter,featurevals,hofx,femp);

options.TolX=acc/10^20;
options.TolFun=acc/10^20;
options.MaxIter=maxiter;
options.MaxFunEvals=maxiter*10;
options.display='none';

[lambda,f,exitflag,output] = minFunc(@(x) FitMaxEntObjFunNew(x,featurevals, hofx, femp),lambda,options);
%keyboard
if output.iterations<=2
%warning('converged suspiciously fast')
options.TolX=acc/10^30;
options.TolFun=acc/10^30;
%options.display='final';
[lambda,f,exitflag,output] = minFunc(@(x) FitMaxEntObjFunNew(x,featurevals, hofx, femp),lambda,options);
end
if output.iterations<=2
%warning('converged suspiciously fast again')
[lambda,f,exitflag,output] = minFunc(@(x) FitMaxEntObjFunNew(x,featurevals, hofx, femp),lambda,options);
%keyboard
end
%nits


[ps,Zfinal]=CalcProbsForFitting(lambda,featurevals,hofx);
[a,b]=MeanVar(ps);
%keyboard
if sum(ps>0)==1 || max(abs(a-mu),abs(b-Sigma))>(1e-5)*n;
%    keyboard
warning('results might be inaccurate');
end

h=lambda(1);
J=lambda(2);
if useBx
   beta=lambda(3);
else
    beta=nan;
end


