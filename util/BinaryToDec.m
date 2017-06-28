function dec = BinaryToDec(bin)
%function dec = BinaryToDec(bin)
%
%convert binary numbers to decimals
%input: a characterstring of 0 and 1s or a logical array
%output: number converted to decimals
%
%example: BinaryToDec('111') gives output 7, BinaryToDec([1 1 1]) likewise

bb = ((size(bin,2)-1):-1:0)';
dec =(bin>0)*(2.^bb);

    
