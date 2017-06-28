function dat=qIsing(dat,H,J,Z,logit)
%compute probability of Ising model with interactions H, couplings J and
%partiion function Z. returns log-probabilities if fourth argument (logit)
%is supplied and true. convention:
% p(x) =1/Z exp(h'x+ 0.5 x'*J*x);


if nargin<=3
    Z=1;
end
if nargin<=4
    logit=0;
end
% unnormalized probability of data under Ising model with H, J
%keyboard


%  keyboard
if max(abs(J))<10^-20;
    dat=dat*H;
elseif size(dat,1)<2^13
    dat=dat*H+.5*diag(dat*J*dat');
else
    % warning('we should add a ridge to get it away from singularity')
    v1= dat*H;
    JJ = sqrtm((J+J')/2);
    %JJ = sqrtm(J);
    A = dat*JJ;
    B = JJ * dat';
    v2=real(sum(A.*conj(B)',2));
    dat=v1+.5*v2;
end

if logit
    dat=dat-log(Z);
else
    dat=exp(dat)/Z;
end
