function [gs,out2]=meancov_2_features(mu, Sigma)
%            [gs] =meancov_2_features(mu, Sigma)
%       [mu,Sigma]=meancov_2_features(gs)
%convert covariance and mean to mean of feature vector, or vice versa

if nargin==2 && nargout==1     % [gs] =meancov_2_features(mu, Sigma)
    %convert covariance and mean to mean of feature vector: Feature vector contains
    %first all mean entries, and then the upper diagonal of all second
    %moments, i.e. covariances with mean product added back in
    Moments=Sigma+mu'*mu;
    gs=[mu(:)',vec(Moments(tril(ones(numel(mu)),-1)==1))'];
    
elseif nargin==1 && nargout==2 % [mu,Sigma]=meancov_2_features(gs)
    %work in opposite direction
    gs=mu;
    dimo=(-1+sqrt(1+8*numel(gs)))/2;
    if rem(dimo,1)>1e-10
        error('number of dimensions does not make sense')
    else
        mu=gs(1:dimo);
        Sigma=zeros(dimo);
        Sigma(tril(ones(dimo),-1)==1)=gs(dimo+1:end);
        Sigma=Sigma+Sigma';
        Sigma=Sigma+diag(mu);
        Sigma=Sigma-mu'*mu;
        gs=mu;
        out2=Sigma;
    end
end

%keyboard