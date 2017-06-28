%simple test script for testing flat models and making sure we get the same
%results in either direction

%first, start by setting up a count model, lets take a DG on 10 states:

mu=0.2;
rho=0.1;
N=10;
[gamma,lambda,count_distrib,model]=fit_flat_dg(mu,rho,10);

states=all_states(10);

ps=count_distrib_2_ps(count_distrib,states);

ps_unsort=count_distrib_2_ps(count_distrib);

count_distrib_2=ps_2_count_distrib(ps,states);

sum(ps)

sum(ps_unsort)

disp('Entropies for 1 2 3 ')

ent(ps)

ent(ps_unsort)

entropy_flat_model(count_distrib)

disp('Log-variances for 1 2 3')

variance_log_probs(ps)

variance_log_probs(ps_unsort)

variance_log_probs_flat_model(count_distrib)

