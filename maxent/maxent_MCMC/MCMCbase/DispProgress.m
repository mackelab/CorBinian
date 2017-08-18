function [s,p,results]=DispProgress(direc)

if nargin==0
    direc=pwd;
else
    cd(direc);
end

s=ReadTrackFile('TrackFile.txt');
load params.mat

N=size(s.metrics,1);

xx=1:N;

Nextrap=max(5*N,5000);

close all,

m=s.metrics;
fl=p.exit_conditions;
m=m./repmat(fl,N,1);
plot(xx,m(1:end,:),'.')


xfit=xx(round(end/4):end);
xfitbig=xx(1):Nextrap;

for k=1:3
    %  keyboard
    poly(k,:)=polyfit((xfit'),log(m(xfit,k)),1);
    trends(k,:)=exp(polyval(poly(k,:),(xfitbig)));
end
hold on
plot((xfitbig),trends)
line([0,Nextrap],[1,1],'color','k')

ab=min(trends<1,[],1);
%keyboard
converged=min(find(ab==1));
if isempty(converged)
    converged=Nextrap;
    titstring='no prediction possible yet, or non-convergence likely';
else
    titstring=['Predicted convergence: ', num2str(converged), ' iterations'];
    if exist('results')==1
        moretime=(now-results.starttime)*24/N*(converged-N);
        titstring=[titstring, ', hours left: ', num2str(moretime)];
    end
end

title(titstring)

xlim([0,converged+10])
ylim([0, ceil(max(m(:)))])
   
set(gca,'ytick',0:ceil(max(m(:))));
if ceil(max(m(:)))<=5
set(gca,'ytick',0:.25:ceil(max(m(:))));
end    
%keyboard





