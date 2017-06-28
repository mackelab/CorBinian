function [P,Lambda]=ProbsDGAnyMarginal(x,supports,gammas,Lambda,acc)


[N,d]=size(x);

[cc, p] = chol(Lambda);
if p > 0
    warning(['Covariance matrix of the latent Gaussian has at least one negative eigenvalue. ' ...
        'Applying Higham-Correction (see help higham).'])
    Lambda = higham(Lambda,1e-10,1e5);
end

gammsleft=0*x;
gammsright=0*x;

for k=1:d
    for kk=1:numel(supports{k})
        gammsright(x(:,k)==supports{k}(kk),k)=gammas{k}(kk);
        if kk>1
            gammsleft(x(:,k)==supports{k}(kk),k)=gammas{k}(kk-1);
        else
            gammsleft(x(:,k)==supports{k}(kk),k)=-inf;
        end
    end
end
%keyboard
P=mvncdf(gammsleft,gammsright,0,Lambda);
