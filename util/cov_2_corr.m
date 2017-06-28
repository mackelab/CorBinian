function A=cov_2_corr(A)
%takes in a covariance matrix,  returns the corresponding correlation
%matrix
%std(A)=diag(A)
%
%JHM 2014

stdA=repmat(sqrt(diag(A)),1,size(A,1));
A=A./stdA./(stdA');