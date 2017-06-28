function A=corr_2_cov_01(A,mu)
%takes in a correlation matrix of binary random variables (0,1) with mean mu,
%returns the corresponding covariance
%matrix
%std(A)=diag(A)
%
%JHM 2014

stdo=sqrt(mu(:).*(1-mu(:)));
stdA=repmat(stdo,1,size(A,1));
A=A.*stdA.*(stdA');