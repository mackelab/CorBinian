function c = countElem(x,sv,ev)

% c = countElem(x,sv,ev)
%  counts the occurences of x for each integer value in [sv,ev), i.e.
%  sv:1:(ev-1)
%
% Code from the paper: 'Generating spike-trains with specified
% correlations', Macke et al., Neural Computation
%
% www.kyb.mpg.de/bethgegroup/code/efficientsampling


if nargin<2
    sv=0;ev=max(x);
end

c =histc(x,sv:ev-1);