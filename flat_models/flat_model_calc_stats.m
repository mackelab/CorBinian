function model=flat_model_calc_stats(count_distrib)


N=numel(count_distrib)-1;
model.N=N;
[model.meancount,model.varcount]=calc_mean_var(count_distrib,0:N);
[model.mean,model.corr]=meanvar_count_2_meancorr(model.meancount,model.varcount,N);
[model.entropy,model.entropy_count]=entropy_flat_model(count_distrib);
[model.var_log_probs,model.mean_log_probs]= flat_model_var_log_probs(count_distrib,exp(1));

model.count_distrib=count_distrib;
