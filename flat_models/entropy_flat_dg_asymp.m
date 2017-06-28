function [H_asymp, H_opper_ubound]    =entropy_flat_dg_asymp(gamma, lambda, Nsubound, Nint,logbase)
%compute asymptotic entropy of flat dichotomized gaussian model with
%gaussian mean gamma and correlation lambda. 

%input: gamma, lambda clear, 
%Nsubound: upper bound calculated from opper.
%Nint: number of integration points:
%logbase: logbase

%output: 
%H_asymp: asymptotic entropy
%H_opper_ubound: opper upper bound (using Nsubound as base
%HC_asump

if nargin<=3 || isempty(Nint)
    Nint=10000;
end

if nargin<=4
    logbase=2;
end
scaler=1/log(logbase);

r=[1:(Nint-1)]/Nint;
%keyboard






H_asymp=-simpint(@(s) integrand(s,gamma,lambda,logbase),-10-gamma,10-gamma, Nint);
%keyboard


if lambda<=10^-10;
    H_asymp=-Phi(gamma).*scaler*log(Phi(gamma))-Phi(-gamma).*scaler*log(Phi(-gamma));
    %H_asymp=H_opper;
%    keyboard
end
if nargout>1
 H_opper_ubound=H_asymp+.5* scaler*log(1+Nsubound*lambda./(1-lambda))./Nsubound;
end
%keyboard


%if nargout>=3
%Pasymp=FlatDGAsymp(gamma, lambda, r);
%Pasymp=Pasymp/sum(Pasymp);
%[HK_asymp,junk,HXgivK_asymp,EntPatterns]=EntropyFromSpikeCount(Pasymp);
%H_asymp=H_asymp/Nint;
%HC_asymp=sum(Pasymp.*(entropy(r).^2))-H_opper^2;
%HC_asymp2=HeatCapacityFromSpikeCount(Pasymp)/Nint^2;
%HC_asymp=simpint(@(s) integrand2(s,gamma,lambda,logbase),-10-gamma,10-gamma, Nint);
%HC_asymp=HC_asymp-H_asymp^2;
%keyboard
%end

function phi=phi(x)
phi=normpdf(x);

function phil=phil(x,lambda)
phil=1/sqrt(lambda)*phi(x/sqrt(lambda));

function Phi=Phi(x)
Phi=normcdf(x);


function invPhi=invPhi(x)
invPhi=norminv(x);

function L=L(s,gamma,lambda)
L=Phi((s+gamma)/sqrt(1-lambda));

function L=log2L(s,gamma,lambda,logbase)
L=normcdfln((s+gamma)/sqrt(1-lambda))/log(logbase);



function p=integrand(s, gamma, lambda,logbase)
LL=L(s,gamma,lambda);
LLL=LL.*log2L(s,gamma,lambda,logbase)+(1-LL).*log2L(-s,-gamma,lambda,logbase);
LLL(isnan(LLL) | isinf(LLL))=0;
p=LLL.*phil(s,lambda);


function p=integrand2(s, gamma, lambda,logbase)
LL=L(s,gamma,lambda);
LLL=LL.*log2L(s,gamma,lambda,logbase)+(1-LL).*log2L(-s,-gamma,lambda,logbase);
LLL(isnan(LLL) | isinf(LLL))=0;
p=(LLL).^2.*phil(s,lambda);

function p=entropy(r,logbase);
scaler=1/log(logbase);
rlogr=r.*log(r)*scaler;
slogs=(1-r).*log(1-r)*scaler;
rlogr(r==0)=0;
slogs(r==1)=0;
p=-rlogr-slogs;

%keyboard



function I=simpint(f,minx,maxx,nsteps);
x=linspace(minx,maxx,nsteps);
stepsize=x(2)-x(1);
y=f(x);
I=sum(y)*stepsize;
%keyboard
