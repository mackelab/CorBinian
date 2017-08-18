/* START FUNCTIONS */
void WriteParams(char *paramsFile,double *lambda,int *first,int *second,int NUMVARS, int NUMFCNS);
void WriteMC(char *mcFile,int **final_mc,int NUMVARS,int NUMDATA);
int WriteMetrics(char *trackFile,double *expec_empir,double *expec_algor,int NUMVARS,int numruns, int tsample, double start_time,double sample_time,double learn_time, double mean_energy_monte_eta, double rel_entropy_monte_theta_lambda, double *exit_conditions);
void expec_diffs(int NUMVARS,double *expec_empir, double *expec_algor,double *expec1_absdiff_algor, double *expec2_absdiff_algor, double *cov_absdiff_algor);
