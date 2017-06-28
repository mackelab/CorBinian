function bin = DecToBinary(dec,d)
%convert decimal numbers to binary ones, output as a logical array.
%
%input: dec the decimal numbers to be converted
%d, the size of the array to output (leave empty to just take minimal array
%size which represents all numbers-- warning: If d is too small leading
%numbers will be truncated and output will be wrong!)

if nargin==2
states=dec2bin(dec,d);
else
    states=dec2bin(dec);
    d=size(states,2);
end

bin=zeros(size(dec,1),d);
for k=1:d
    bin(:,k)=str2num(states(:,k));
end
bin=logical(bin);
    
