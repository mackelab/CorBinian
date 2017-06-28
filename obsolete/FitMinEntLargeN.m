function [pmf,fval,exitflag,output,lambda]=FitMinEntLargeN(N,meano, varo, flipit)

%x = linprog(f,A,b,Aeq,beq,lb,ub,x0,options)
%maximize f' x subject to
% Ax leq b
% Aeq x=beg
% lb leq x leq ub

%equality matrices:
myconst=ones(N+1,1);
myx=(0:N)';
myxsquare=((0:N).^2)';

Aeq=[myconst, myx/N, myxsquare/N^2]';
%beq=[1,meano*N, meano*(1-meano)*N+varo*N*(N-1)+meano^2*N^2]';
beq=[1,meano, varo+meano^2];
ops=optimset('display','off');
lb=zeros(N+1,1);
lognchoosek=gammaln(N+1)-gammaln((0:N)'+1)-gammaln(N-(0:N)'+1);
f=lognchoosek*log2(exp(1));
if nargin==4 && flipit
    f=-f;
end
%keyboard
%keyboard
[pmf,fval,exitflag,output,lambda]=linprog(f,[],[],Aeq, beq, lb, [],[],ops);

lambda.Aeq=Aeq;
lambda.beq=beq;
lambda.c=f;
lambda.x=pmf;
lambda.fval=fval;
