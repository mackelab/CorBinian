function [PsMarg,PsInd]=CalcMarginals(states,Ps)
%function [PsMarg,PsInd]=CalcMarginals(states,Ps)
%
%given binary vectors 'states' (as logical array, one line for each state)
%and the probability mass of each state in the vector Ps, calculates the
%marginal distribution over each binary variable (i.e. the probability that
%this neuron fires or that this spin is positive. i.e. P(x_i)= 1 for each
%i.
%Optinally also returns the probability of each state under the assumption
%that the different dimensions are independent.
%


if min(size(states))==1
    %if states are given as decimals, convert to binary first:
    states=DecToBinary(states);
end

d=size(states,2); %dimensionality of states
N=size(states,1); %number of states

if nargin==1 %assume uniform distribution if Ps is missing
    Ps=ones(N,1)/N;
end


for k=1:d
    PsMarg(k)=sum(Ps(states(:,k)==1));
end

if nargout==2
    PsInd=CalcIndep(PsMarg,states);
end