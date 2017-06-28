#ifndef LEARN_H
#define LEARN_H



/* START GLOBAL VARIABLES */
/* start iteration parameters */
int NUMVARS;
int NUMFCNS;
int NUMDATA;
int NUMSAMPS;
/* done iteration parameters */

/* start Monte Carlo samples */
int *count_monte;
int *bound_monte;
int **sInds_monte;
/* lambda */
double *lambda;
double *g_monte_lambda;
double Z_monte_lambda;
/* eta */
double *eta;
double *g_monte_eta;
double Z_monte_eta;
double *sum_fjexpg_monte_eta;
/* done Monte Carlo samples */

/* start positive examples */
int *count_pos;
int *bound_pos;
int **sInds_pos;
/* theta */
double *g_pos_theta;
/* done positive examples */

/* start checks */
double rel_entropy_monte_theta_lambda;
    double sum_gexpg_monte_eta;
double mean_energy_monte_eta; /* energy = g_monte_lambda */
    double sum_energyexpg_monte_eta;
double loss_pos;
    double sum_g_pos_theta;
double old_loss_pos;
/* done checks */

/* start expectations */
int **sigma;
double expec_empir[50000];
double init_lambda[50000];
/* done expectations */

/* start counters */
int *first;
int *second;
int **combine;
/* done counters */
/* DONE GLOBAL VARIABLES */

#endif
