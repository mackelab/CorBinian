function [abc,avec,bmat,cmat]=ApproxFinishLines(N,mu,Sig)
%calcualte approximate standard deviations for the empirical mean, cov, and correlation o
%of binary data with mean mu and covariance Sig, assuming a datasets of
%size N. N should either be the size of the MC sample, or the size of the
%input dataset. In principle though, the MC sample should always be bigger
%than the original dataset.
%calculates it for every dimension separately, and then averages across
%dimensions, but also returns the separate entries in avec bmat, cmat.
%keyboard
d=numel(mu);

mu=mu(:)';

avec=(mu.*(1-mu))/N;
abc(1)=sqrt(mean(avec));

Moment=Sig+repmat(mu,d,1).*repmat(mu,d,1)';
bmat=Moment.*(1-Moment)/N;
abc(2)=sqrt(mean(vec(bmat)));

cmat=Moment./(repmat(mu,d,1).*repmat(mu,d,1)')/N;

abc(3)=sqrt(mean(vec(cmat)));
