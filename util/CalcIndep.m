function [PsInd]=CalcIndep(mu,states)
%function [PsInd]=CalcIndep(mu,states)
%
%calc independent distribution over binary patterns assuming means mu
%
mu=mu(:);

PsInd=exp(states*log(mu)+(1-states)*log(1-mu));