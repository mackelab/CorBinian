function [HC,EntModel,EntPatterns,HModelSq]=HeatCapacityFromSpikeCount(ps,logbase)
%calculates the entropy of a spike count histogram "ps" and the entropy of
%the underlying population model, under the assumption that all patterns
%with the same number of spike have the same probability. This is the
%maximum entropy model consistent with a given spike count histogram. For
%Ising models, it corresponds to an Ising model with constant pairwise
%couplings. 
%keyboard

%NOTE: THIS FUNCTION MIGHT HAVE A SMALL BUG HAVING TO DO WITH LOG2. I used
%log2 for the result in the PRL paper, but this might be wrong
if nargin==1 
    logbase=2;
end
scaler=1/log(logbase);
    
    
N=numel(ps)-1;
ps=ps(:);
plogp=ps.*(log(ps)*scaler);
plogp(ps==0)=0;

EntCount=-sum(plogp);
lognchoosek=gammaln(N+1)-gammaln((0:N)'+1)-gammaln(N-(0:N)'+1);


lognchoosek=lognchoosek*scaler;


EntPatternsgivenModel=sum(ps.*lognchoosek);
EntModel=EntCount+EntPatternsgivenModel;
EntPatterns=lognchoosek;
%keyboard

sqlogdif=(log(ps)*scaler-lognchoosek).^2;
psqlogdif=ps.*sqlogdif;
psqlogdif(ps==0)=0;
psqlogdif(isnan(psqlogdif))=0;


HModelSq=sum(psqlogdif);

HC=HModelSq-EntModel^2;
%keyboard
