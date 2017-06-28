% This script generates data means E[f(X)] from P(K) under the assumption
% of the flat model. 
% It requires a (N_max+1)-by-M matrix called 'pcounts' containing the
% P(K)s for M different subpopulation sizes with the largest subpopulation
% size being N_max. 
% The script compute_Heat_Curves_tidy.m procudes such a matrix.

Efx = cell(size(pcounts,2),1);

figure; 
clrs = jet(size(pcounts,2)); 
for i = 1:length(Efx);

% get mean firing rate, mean covariance, mean E[xi * xj]
model=flat_model_calc_stats(pcounts(1:Ns(i)+1, i));
[~,~,cov]=meanvar_count_2_meancorr(model.meancount,model.varcount,model.N);
E1 = model.mean;
E2 = cov + E1^2;

Efx{i} = zeros( model.N * (model.N  + 3) / 2 + 1 , 1 );
Efx{i}(1:model.N) = E1;
Efx{i}(model.N + (1:model.N*(model.N-1)/2)) = E2;
Efx{i}(model.N*(model.N+1)/2+1 : end) = pcounts(1:Ns(i)+1, i);
plot(Efx{i},'color', clrs(i,:)); hold on
end

