function S = DecToBinary_alt(X,n)
% converts numbers in row vector into a binary string of length n
% alternative to DecToBinary.m-- not sure which one is better, have not
% checked yet. Both give different results if n is too small!! This
% function is not used anywhere in the code, just added it in in case we
% want to go back to it at some point.


cc = zeros(size(X));
S = zeros(n,size(X,2));
for i=n:-1:1
    idx = X>=2^(i-1);
    S(n-i+1,idx)=1;
    cc(idx) = cc(idx)+1;
    X(idx) = X(idx) - 2^(i-1);
end

S=S';