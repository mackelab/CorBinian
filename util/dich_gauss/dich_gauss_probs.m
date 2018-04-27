function [P]=ProbsDG(x,gamma,Lambda,varargin)
[N,d]=size(x);

%tic

gammsleft=zeros(N,d);
gammsright=zeros(N,d);

gammsleft(x<=0)=-inf;
gammsright(x>0)=+inf;

chunksize=1000;
numchunks=ceil(N/chunksize);

gamma=gamma(:)';

for k=1:numchunks

indie=(k-1)* chunksize+1:min(N,k*chunksize);
P(indie)=mvncdf(gammsleft(indie,:),gammsright(indie,:),gamma,Lambda,varargin);
end

%toc

