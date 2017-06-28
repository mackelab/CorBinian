function fx = maxEnt_features_slim(x, model)
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
    case 'ising_count_l_0' % same as 'ising_count', but has feature for K=0
        pairs=nchoosek(1:d,2);
        count_indicators=zeros(d+1,1);
        count_indicators(sum(x)+1) = 1; 
        %description=[[1:d; (1:d)*nan],pairs,[0:d; (0:d)*nan]];  
        fx=[x;x(pairs(:,1)).*x(pairs(:,2));count_indicators];
end