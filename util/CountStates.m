function [freqs,PNoHit,states]=CountStates(samples,states,P)
%count frequency (i.e. histogram) over samples. If used with one argument,
%just calculates histogram over all unique states, if 'states' is supplied
%as an argument, only considers samples which match one of the state. If
%third argument P is also supplied, re-weights histogram by the a-priori
%mass of each state.

%input: samples: an array of binary vectors (one line of each sample) or a vector of decimal numbers, the frequency of which we want to count 
%       states (optional): an array of binary vectors (or a vector of decimals) of all the possible states to be included. (Set empty if all unique samples should be taken)  
%       P (optional): The a priori probability mass function of each sample. Leave
%       empty to take uniform. 
%
%output: 
%freqs: A vector with as many entries as 'states'. For each state, tells us
%how often it occured in the samples. 
%PNoHit: The total probability mass/frequency of samples which did not match any
%of the states
%states: The list of all possible states (if states is used as input, this
%is just the input passed out again)

if islogical(samples)
    samples=BinaryToDec(samples);
end

if nargin==1 || isempty(states)
    states=unique(samples);
end

if islogical(states)
    states=BinaryToDec(states);
end

freqs=zeros(size(states,1),1);
N=numel(samples);

if nargin<=2
    for k=1:numel(states)
        freqs(k)=nnz(samples==states(k));
    end
    freqs=freqs/N;
else
    for k=1:numel(states)
        freqs(k)=sum(P(samples==states(k)));
    end  
    freqs=freqs/sum(P);
end

PNoHit=1-sum(freqs);



states=DecToBinary(states);
