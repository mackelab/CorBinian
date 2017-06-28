function A=AllStates(d);
%make list of all binary states with d elements ordered in ascending order

list=[0:(2^d-1)]';

A=DecToBinary(list);
