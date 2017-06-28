function [Z,states,allstates,mu,CC,Ps]=StupidCompZ(H,J);


d=numel(H);
Nstates=2^d;
states=dec2bin((0:Nstates-1)');
allstates=zeros(Nstates,d);
for k=1:d
    allstates(:,k)=str2num(states(:,k));
end
%keyboard
Ps=qIsing(allstates,H,J);
%keyboard
Z=sum(Ps);
Ps=Ps/Z;

if nargout>3
    mu=sum(repmat(Ps,1,d).*allstates);
    mystates=allstates-repmat(mu,Nstates,1);
    CC=(repmat(Ps,1,d).*mystates)'*mystates;
end

allstates=logical(allstates);
