function s=DecToBin(d,n)
%DECToBIN Convert decimal integer to a binary vector.
%   DEC2BIN(D) returns the binary representation of D as a string.
%   D must be a non-negative integer smaller than 2^52.
%second argument n uses the minimum number of digits to use
%

if nargin==1
    n=1;
end

[f,e]=log2(max(d)); %#ok How many digits do we need to represent the numbers?
s=(rem(floor(d*pow2(1-max(n,e):0)),2));
