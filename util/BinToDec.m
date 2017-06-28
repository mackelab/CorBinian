function d=BinToDec(b)
%convert binary vector {0,1}s to decimal numbers
%terrible function using matlab-builtin functions that are based on
%strings, could probably be sped up by a lot if needed

b=b==1;
bb=repmat('0',size(b));
bb(b)='1';

d=bin2dec(bb);
