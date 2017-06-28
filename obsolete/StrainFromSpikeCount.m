function [strainlog2,strainloge]=StrainFromSpikeCount(p);

%calculate strain from spike count distribution of three states, really
%hacky funciton, but me does not care

%patterns are 
if numel(p)~=4
    error('asdfsdf')
end

P000=p(1);
P001=p(2)/3;
P010=p(2)/3;
P011=p(3)/3;
P100=p(2)/3;
P101=p(3)/3;
P110=p(3)/3;
P111=p(4);

strainlog2=1/8*log2(P100*P010*P001*P111/(P000*P011*P101*P110));
strainloge=1/8*log(P100*P010*P001*P111/(P000*P011*P101*P110));



%strain is log(p100 p010
