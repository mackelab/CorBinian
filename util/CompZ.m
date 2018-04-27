function [Z,p,x,m,C]=CompZ(h,J,x)

% [Z,p,x,m,C] = CompZ(h,J)
%   Computes the normalized probabilities p under the Ising model with
%   parameters h and J by brute force enumeration of all possible states.
%   Implementation avoids two computational bottlenecks: 1) We generate the
%   2^D binary states in O(D) instead of O(2^D) and 2) we avoid memory
%   problems when computing the probablities by factorizing the involved
%   matrizes into managable junks.
%
%   p       probabilities
%   x       binary states
%   m       mean vector
%   C       covariance matrix
%   Z       normalizing constant
%
%   h       activation vector
%   J       interaction matrix
%
% PHB 2007-11-01
% modified JHM 2014

d = length(h);

% compute all possible binary patterns


if nargin==3
    %
else
   % keyboard
    x = DecToBinary((0:2^d-1)',d);
end

% compute Ising probabilites
p = qIsing(x,h,J);

% normalize
Z=sum(p);
p=p/Z;

% compute mean and covariance, if asked for
%keyboard
if nargout>3
    chunksize=5000;
    numchunks=ceil(numel(p)/chunksize);
    if numchunks==1
        m=sum(repmat(p,1,d).*x);
        xx=x-repmat(m,2^d,1);
        C=(repmat(p,1,d).*xx)'*xx;
    else
        m=0;
        C=0;
        for k=1:numchunks
            indie=chunksize*(k-1)+1:min(chunksize*k,numel(p));
            locx=x(indie,:);
            m=m+sum(repmat(p(indie),1,d).*locx);
            xx=locx-repmat(m,numel(indie),1);
            C=C+(repmat(p(indie),1,d).*xx)'*xx;
        end
    end
end

x=logical(x);






