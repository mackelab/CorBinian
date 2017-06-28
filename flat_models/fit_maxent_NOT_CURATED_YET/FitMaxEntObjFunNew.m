function    [L,dL]=FitMaxEntObjFunNew(lambda,featurevals,hofx,femp)

s=lambda'*featurevals+hofx;

Q=exp(s);
Z=sum(Q);
logZ=log(Z);
P=Q/Z;

L=logZ-lambda'*femp;

%keyboard
dL= featurevals*P'-femp;

