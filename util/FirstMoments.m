function [moments]=FirstMoments(ps,r,m)
%calculate first m moments of distribution ps on support r

ps=ps/sum(ps);

if nargin==1
    n=numel(ps)-1;
    r=[0:n]';
end

if nargin<=2
    m=4;
end


for k=1:m
    moments(k)=sum(r.^k.*ps(:));
end
