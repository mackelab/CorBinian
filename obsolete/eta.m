function eta=eta(r,logbase)

if nargin==1
    logbase=2;
end

if logbase==2
eta=-r.*log2(r)-(1-r).*log2(1-r);
else
eta=-r.*log(r)-(1-r).*log(1-r);
end    

eta(isnan(eta))=0;
eta(eta>=inf)=0;
