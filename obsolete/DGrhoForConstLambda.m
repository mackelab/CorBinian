function rhotrans=DGrhoFOrConstLambda(r, lambda)

if nargin==1
    lambda=1/2;
end

vars=r.*(1-r);
gammas=norminv(r);
%keyboard
rhotrans=(mvncdf([gammas(:),gammas(:)],[0,0], [1,lambda;lambda,1])-r(:).^2)./vars(:);
rhotrans(isnan(rhotrans))=0;
