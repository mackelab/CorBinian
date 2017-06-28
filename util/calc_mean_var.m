function [M,V]=calc_mean_var(ps,r)
%given a probability mass function ps, calculate its mean M and variance V.
%If second argument r is supplied, r gives the possible states of the
%random variable. If r is not supplied, r=[0:n] is assumed, if r==1, then
%r=[0:n]/n$ is assumed.
%
%JHM 2014

ps=ps/sum(ps(:));

if min(size(ps))==1 || min(size(r))==1
    if nargin==1
        n=numel(ps)-1;
        r=[0:n]';
    elseif nargin==2 && numel(r)==1 && r==1;
        n=numel(ps)-1;
        r=[0:n]'/n;
        
    end
   % keyboard
    M=sum(r(:).*ps(:));
    V=sum(ps(:).*(r(:)-M).^2);
elseif min(size(ps))>=2 || min(size(r))>=2
    warning('Undocumented and untested behaviour in 2-dimensional case')
%    keyboard
       n=size(ps)-1;
       if nargin==1
       [r1,r2]=meshgrid(0:n(1),0:n(2));
       elseif nargin==2 && numel(r)==1
       [r1,r2]=meshgrid((0:n(1))/n(1),(0:n(2))/n(2));
       else
           r1=r(:,1);
           r2=r(:,2);
       end
       
       M(1)=sum(ps(:).*r1(:));
       M(2)=sum(ps(:).*r2(:));
       V(1,1)=sum(ps(:).*r1(:).*r1(:))-M(1)^2;
       V(1,2)=sum(ps(:).*r1(:).*r2(:))-M(1)*M(2);
       V(2,2)=sum(ps(:).*r2(:).*r2(:))-M(2)^2;
       V(2,1)=V(1,2);

end
