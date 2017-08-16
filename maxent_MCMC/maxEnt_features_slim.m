function fx = maxEnt_features_slim(x, model)
% calculature features for max-ent model fitting


% first argument: the vector to which the feature maps should be supplied
% inputs: 
%  x: a binary vector of length d
% map: either integers 1, 2, 3, ('first', 'second', 'third-order' maxent),
% or 'ising' (same as 2) or 'ising_count' for ising model with additional 
% constraints on the total activity count
% or a function handle for using user-supplied feature maps
%
%outputs:
% fvals, the calculated features. If map=1, then same as x, if map=2, this
% is a copy of x concatenated with all two-tupels without repetion, if
% map=3, then this is further concatenated with all three-tupels without
% any repetions. (Note that this is NOT what you might want if x is
% different from a binary (0,1) repetiation, as it excludes all the x.^2,
% and x.^3 features etc, as they are redundant in this representation!

% function adapted from the pop_spike code base 
% https://bitbucket.org/mackelab/pop_spike

 d = length(x);
 switch model
    case {2,'ising'}
        pairs=nchoosek(1:d,2);
        %description=[[1:d; (1:d)*nan],pairs];  
        fx=[x;x(pairs(:,1)).*x(pairs(:,2))];
    case 3
        pairs=nchoosek(1:d,2);
        triplets=nchoosek(1:d,3);
        %description=[[1:d; (1:d)*nan;(1:d)*nan],[pairs; pairs(1,:)*nan],triplets];  
        fx=[x;x(pairs(:,1)).*x(pairs(:,2));x(pairs(:,1)).*x(pairs(:,2)).*x(triplets(:,3))];
    case 'ising_count'
        pairs=nchoosek(1:d,2);
        count_indicators=zeros(d,1);
        sum_x = sum(x);
        if sum_x>1
         count_indicators(sum_x)=1;
        end
        %description=[[1:d; (1:d)*nan],pairs,[1:d; (1:d)*nan]];  
        fx=[x;x(pairs(:,1)).*x(pairs(:,2));count_indicators];
        %keyboard
    case 'k_pairwise' % same as 'ising_count', but has feature for K=0
        pairs=nchoosek(1:d,2);
        count_indicators=zeros(d+1,1);
        count_indicators(sum(x)+1) = 1; 
        %description=[[1:d; (1:d)*nan],pairs,[0:d; (0:d)*nan]];  
        fx=[x;x(pairs(:,1)).*x(pairs(:,2));count_indicators];
end