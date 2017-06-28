function hc = binHist(S)
% hc = binHist(S)
%   Generates histogram for binary vectors in S
%   S = N * D matrix
%
% Code from the paper: 'Generating spike-trains with specified
% correlations', Macke et al., Neural Computation
%

S=S';
D = size(S,1);
C = binBinaryToDec(S);
hc = countElem(C,0,2^D)';
