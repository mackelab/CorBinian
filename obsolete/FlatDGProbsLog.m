function P=FlatDGProbs(gamma,lambda,k,n,intsteps, logit)


%keyboard
if nargin==4
    intsteps=5000;
end
if nargin<=5
    logit=0;
end




for i=1:numel(k)

    r=k(i)/n;
    if k(i)>0 & abs(lambda)>.00001 & n>50 & k(i)<n;
        %from the asymptotic theory, we know that the integral is concentrated
        %around
        so=sqrt(1-lambda) .* invPhi(r)-gamma;
        %and all the relevant mass in on a scale of
        %keyboard
        intrange=1/(phi(invPhi(r))/sqrt(1-lambda)*sqrt(n+500)/sqrt(500)/sqrt(r*(1-r)));
        keyboard
    else
        so=-gamma;
        intrange=5;
    end
    %keyboard
    %plot(linspace(so-intrange, so+intrange,intsteps),integrand(linspace(so-intrange, so+intrange,intsteps),r,n,gamma,lambda))
    %keyboard
    p(i)=(logsimpint(@(s) integrand(s,r,n,gamma,lambda),so-intrange,so+intrange, intsteps));
    %s=linspace(so-intrange,so+intrange,intsteps);
    %integrand(s,r,n,gamma,lambda)
    %keyboard
    %pp(i)=log(quad(@(s) integrand(s,r(i),n,gamma,lambda),-5,5));
    noverk(i)=-betaln(n-k(i)+1,k(i)+1)-log(n+1);
    P(i)=p(i)+noverk(i);

    %    keyboard
end
if ~logit
    P=exp(P);
end

%
%keyboard

%keyboard

function phi=logphi(x)
phi=normpdfln(x);

function phi=phi(x)
phi=normpdf(x);


function phil=logphil(x,lambda)
phil=-0.5*log(lambda)+logphi(x/sqrt(lambda));


function Phi=Phi(x)
Phi=normcdf(x);

function LogPhi=LogPhi(x);
LogPhi=normcdfln(x);



function invPhi=invPhi(x)
invPhi=norminv(x);


function L=LogL(s,gamma,lambda)
L=LogPhi((s+gamma)/sqrt(1-lambda));



function p=integrand(s, r, n, gamma, lambda)
LL=LogL(s,gamma,lambda);
LLminus=LogL(-s,-gamma,lambda);
LLL=(r*n).*LL+(n-r.*n)*LLminus;
diffo=log(s(2)-s(1));
LLL(isnan(LLL))=0;
p=logphil(s,lambda);
checksum=logsumexp(p')+diffo;
p=p-(checksum);
checksum=logsumexp(p')+diffo;
p=p-(checksum);
%keyboard
if abs(1-exp(logsumexp(p')+log(s(2)-s(1))))>1e-15;
    checksum=logsumexp(p')+diffo;
    p=p-(checksum);
    if abs(1-exp(logsumexp(p')+log(s(2)-s(1))))>1e-15;
        checksum=logsumexp(p')+diffo;
        p=p-(checksum);

        if abs(1-exp(logsumexp(p')+log(s(2)-s(1))))>1e-15;
            checksum=logsumexp(p')+diffo;
            p=p-(checksum);
            warning('integration in DG might be inaccurate');
            %keyboard
        end
    end
end
p=p+LLL;


function I=logsimpint(f,minx,maxx,nsteps);
x=linspace(minx,maxx,nsteps);
stepsize=x(2)-x(1);
y=f(x);
%keyboard
I=logsumexp(y')+log(stepsize);

%keyboard






