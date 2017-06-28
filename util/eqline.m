function [h]=eqline(a);
%print line of equality into a plot
%
% JHM 2014

if nargin==0
    a=gca;
end

[x]=xlim(a);
[y]=ylim(a);
mini=min(x(1),y(1));
maxi=max(x(2),y(2));

hh=line([mini,maxi],[mini,maxi],'color','k');

if nargout==1
    h=hh;
end
