function [hists,ranges,joints2D]=EstimateDiscreteMarginal(A,includezeros);
%a simple function for estimating the marginals of states in a
%matrix A which is of size num trials by num variables. The range of each dimension is taken to be as large
%as needed.
%very slow and inefficient, but I dont care (yet)

[N,d]=size(A);
for k=1:d
    %    keyboard
    if nargin==1 || includezeros==false
        [ranges{k}]=unique(A(:,k));
    else
        ranges{k}=(0:max(A(:,k)));
    end
    rangesize(k)=numel(ranges{k});

end
J=num2cell(j);
numstates=prod(rangesize);


for k=1:d
    for kk=1:numel(ranges{k});
        hists{k}(kk)=sum(ranges{k}(kk)==A(:,k));
    end
    hists{k}=hists{k}/N;
    if nargout>1
        for kk=k+1:d
            for kkk=1:numel(ranges{k})
                for kkkk=1:numel(ranges{kk})
                    joints2D{k,kk}(kkk,kkkk)=sum(ranges{k}(kkk)==A(:,k) & ranges{kk}(kkkk)==A(:,kk));
                end
            end
            joints2D{kk,k}=joints2D{k,kk}/sum(vec(joints2D{k,kk}));
            joints2D{kk,k}=joints2D{k,kk}';
        end
    end
end


%keyboard
