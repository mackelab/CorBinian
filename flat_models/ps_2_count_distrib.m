function [pK]=ps_2_count_distrib(states, ps)
%convert full probability mass function over states to count distribuion.
%Counts are defined as K = sum(states,2) and for binary states range
%between 0 and d. 
%
% Input: 
%  - states: N-by-d matrix of N states in d-dimensional space. 
%  - ps:     N-dim. vector of probabilities for each state
% Output:
%  - pK: (N+1)-dim. vector of probability of counts K = 0,...,d
%
%also see CountOnes

if ~islogical(states) && any(any(states>1)); %
    states=DecToBinary(states);              % ? ... just don't touch it...
end                                          %

K = sum(states,2);     % compute counts
Kmax = size(states,2); % maximum possible count

pK = zeros(1, Kmax+1);
if nargin==1 % use empirical distribution for ps
    for k=0:Kmax
        pK(k+1)=sum(K==k);
    end
    
else         % make use of provided ps
    for k=0:Kmax
       pK(k+1)=sum(ps(K==k));
    end
end

pK=pK/sum(pK); % normalize (should not be needed if nargin > 1...)
