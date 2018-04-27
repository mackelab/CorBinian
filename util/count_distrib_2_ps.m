function [ps]=count_distrib_2_ps(p_count,states)
%convert a spike-count distribution p_count to a probability distribution
%over all states, assuming a 'flat' model, i.e. that all states with the
%same number of spikes are equally likely. This functon can work in 2
%modes: If second argument 'states' is supplied, then it sorts ps according
%to the states in the matrix states, i.e. such that ps(k)= P(states(k,;)).
%If there is no second argument, it just lists them in order of increasing
%k.


N=numel(p_count)-1;

if nargin==1 %unsorted mode, just make long list of ps, sorted by k:
    lognchoosek=gammaln(N+1)-gammaln((0:N)'+1)-gammaln(N-(0:N)'+1);
    mynchoosek=round(exp(lognchoosek));
    ps=[];
    for k=1:(N+1);
        ps_loc=p_count(k)/(mynchoosek(k));
        ps=[ps; repmat(ps_loc,(mynchoosek(k)),1);];
    end
    %keyboard
elseif nargin==2 %sorted mode, sort elements of ps according to the list of states

   state_counts=sum(states,2);
   ps=zeros(size(state_counts));
   for k=1:N+1
       index=(state_counts==(k-1));
       ps(index)=p_count(k)/sum(index);
   end
end