#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "expec_err.h"

void WriteParams(char *paramsFile,double *lambda,int *first,int *second,int NUMVARS, int NUMFCNS) {  
	FILE *fparams;
	int j;
	char paramsFileLoc[300];
	strcpy(paramsFileLoc,paramsFile);
	strcat(paramsFileLoc,".txt");

	if((fparams = fopen(paramsFileLoc, "w")) == NULL) {
        fprintf(stderr, "Bad file name: %s\n", paramsFile);
        exit(-1);
    }
    fprintf(fparams, "%d\n", NUMVARS);
    for(j = 0; j < NUMVARS; j++) {
        fprintf(fparams, "%d %.14f\n", j, lambda[j]);
    }
    for(j = NUMVARS; j < NUMFCNS; j++) {
        fprintf(fparams, "%d %d %.14f\n", first[j], second[j], lambda[j]);
    }
	fclose(fparams);
}

void WriteMC(char *mcFile,int **final_mc,int NUMVARS,int NUMDATA) {
	FILE *fmc;
	int j,mc;
    if((fmc = fopen(mcFile, "w")) == NULL) {
        fprintf(stderr, "Bad file name: %s\n", mcFile);
        exit(-1);
    }
    fprintf(fmc, "%d\n", NUMVARS);
    for(j = 0; j < NUMVARS; j++) {
        for(mc=0; mc < NUMDATA;mc++) {
            fprintf(fmc, "%d\n",final_mc[j][mc]);
        }
    }
	fclose(fmc);
}

int WriteMetrics(char *trackFile,double *expec_empir,double *expec_algor,int NUMVARS,int numruns, int tsample, double start_time, double sample_time, double learn_time, double mean_energy_monte_eta, double rel_entropy_monte_theta_lambda, double *exit_conditions) {
	FILE *ftrack;
	double expec1_absdiff_algor, expec2_absdiff_algor, cov_absdiff_algor;
	double elapsed_time;
	double convergence_gap[4];
	int success[4];
	int stop_me;
	
	if((ftrack = fopen(trackFile, "a")) == NULL) {
        fprintf(stderr, "Bad file name: %s\n", trackFile);
        exit(-1);
    }
	
	expec_diffs(NUMVARS, expec_empir, expec_algor,&expec1_absdiff_algor, &expec2_absdiff_algor, &cov_absdiff_algor);
	
	elapsed_time = (double) clock() - start_time;
	
	convergence_gap[0]=expec1_absdiff_algor-exit_conditions[0];
	convergence_gap[1]=expec2_absdiff_algor-exit_conditions[1];
	convergence_gap[2]=cov_absdiff_algor-exit_conditions[2];
	convergence_gap[3]=exit_conditions[3]-elapsed_time;
	
	
	success[0]=convergence_gap[0]<0;
	success[1]=convergence_gap[1]<0;
	success[2]=convergence_gap[2]<0;
	success[3]=convergence_gap[3]<0;
	/*fprintf(stderr,"%f %f %f %d \n", expec1_absdiff_algor, exit_conditions[0], convergence_gap[0], success[0]);*/
	
	stop_me=(((success[0]==1) &&  (success[1]==1) && (success[2]==1)) || success[3]==1);


	
	/*mean_energy_monte_eta=1234567;*/
	/*fprintf(ftrack, " %.14f %.14f %.14f %d %d %.14f %.14f %.14f \n", expec1_absdiff_algor, expec2_absdiff_algor, cov_absdiff_algor,numruns, tsample, ((double)elapsed_time)/((double) CLOCKS_PER_SEC),rel_entropy_monte_theta_lambda,mean_energy_monte_eta);*/
	
	fprintf(ftrack, "%d %d %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %d %d %d %d \n", numruns, tsample, expec1_absdiff_algor, expec2_absdiff_algor, cov_absdiff_algor, ((double)elapsed_time)/((double) CLOCKS_PER_SEC),((double)sample_time)/((double) CLOCKS_PER_SEC),((double) learn_time)/((double) CLOCKS_PER_SEC),mean_energy_monte_eta,rel_entropy_monte_theta_lambda, success[0], success[1], success[2], success[3]);
	fclose(ftrack);
	
	return stop_me; 
}

