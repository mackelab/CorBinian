/*
 * Uses learning rate to change coeffs at each time step.
 * lambda + eta = theta (new lambda)
 * 
 * COMPILE: make all
 *  */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "sampler.h"
#include "expec_err.h"
#include "learn_iterative.h"
#include "learn_stop.h"
#include "OutputFunctions.h"
#include "ReadInputs.h"

/*#define DEBUG*/

/* number of samples generated to burn-in the Monte Carlo Markov Chain */
#define BURNINLENGTH 0
/* random number generator seed */
#define SEED 2


int main(int argc, char *argv[]) {
    /* START DEFINING */
    /* start iteration parameters; */
    int TCOUNT;
    int NUMRUNS;
    /* done iteration parameters */
    
    /* start final monte carlo sample */
    int **final_mc;
    
    /* start checks */
    double *lambda_backup;
    /* done checks */
    
    /* start expectations */
    double *expec_algor;
	
    double Z_algor_frac_accum_algor_immed;
    /* done expectations */
        
    /* start files */
    char trackFile[300];
    char paramsFile[300];
    char mcFile[300];
	char numrunsString[30];
	char paramsFileTemp[300]; 
    
    /* start counters, miscellaneous */
	/* clock watchers */
    double start_time, learn_time, sample_time;
    int mc, j, i0, i1;
    int numruns, tsample;
	double exit_conditions[4];
	int stop_me;
    /* recalculate time 
    int recalc_expec;
	recalc_expec=0;*/
    
    /*some variables for timing the duration of each learning and sampling step*/
    learn_time=0;
    sample_time=0;
    
	/* DONE DEFINING */

    /* START INITIALIZING */
	
	/*****************************************************************************/
    /*************************** START READING INPUTS  ***************************/
	/*****************************************************************************/
	
   Inputs(argv,argc, &NUMVARS, &NUMRUNS, &TCOUNT, &NUMSAMPS, &NUMDATA,trackFile,paramsFile,mcFile,expec_empir,init_lambda,exit_conditions);
	
    /***************************** END READING INPUTS  ***********************/
    
    /*****************************************************************************/
    /*************************** START ALLOCATING SPACE  *************************/
	/*****************************************************************************/
	
    NUMFCNS = NUMVARS + NUMVARS * (NUMVARS - 1) / 2;
    
    /* start Monte Carlo samples */
    count_monte = (int *) malloc(NUMFCNS * sizeof(int));
    bound_monte = (int *) malloc(NUMFCNS * sizeof(int));
    sInds_monte = (int **) malloc(NUMFCNS * sizeof(int *));
    for(j = 0; j < NUMFCNS; j++) {
        bound_monte[j] = 1;
        sInds_monte[j] = (int *) malloc(bound_monte[j] * sizeof(int));
        assert(sInds_monte[j] != NULL);
    }
    /* lambda */
    lambda = (double *) malloc(NUMFCNS * sizeof(double));
    g_monte_lambda = (double *) malloc(NUMSAMPS * sizeof(double));
    /* eta */
    eta = (double *) malloc(NUMFCNS * sizeof(double));
    g_monte_eta = (double *) malloc(NUMSAMPS * sizeof(double));
    sum_fjexpg_monte_eta = (double *) malloc(NUMFCNS * sizeof(double));
    /* done Monte Carlo samples */
    
    /* start final Monte Carlo sample */
    final_mc = (int **) malloc(NUMVARS * sizeof(int *));
    for(i0 = 0; i0 < NUMVARS; i0++) {
        final_mc[i0] = (int *) malloc(NUMDATA * sizeof(int));
    }
    /* done final Monte Carlo sample */
     
    
    /* start expectations */
    sigma = (int **) malloc(NUMVARS * sizeof(int *));
    for(i0 = 0; i0 < NUMVARS; i0++) {
        sigma[i0] = (int *) malloc(NUMDATA * sizeof(int));
    }
    /* done expectations */

    /* start checks */
    lambda_backup = (double *) malloc(NUMFCNS * sizeof(double));
    expec_algor = (double *) malloc(NUMFCNS * sizeof(double));
    /* done checks */
    
    /* start counters */
    first = (int *) malloc(NUMFCNS * sizeof(int));
    second = (int *) malloc(NUMFCNS * sizeof(int));
    combine = (int **) malloc(NUMVARS * sizeof(int *));
    for(i0 = 0; i0 < NUMVARS; i0++) {
        combine[i0] = (int *) malloc(NUMVARS * sizeof(int));
    }
    /* done counters */
	
    /**********************done allocating space *****************/
    
    /********************** start setting counters **********************/
    
	j = NUMVARS;
    for(i0 = 0; i0 < NUMVARS; i0++) {
        combine[i0][i0] = -1;
        for(i1 = i0 + 1; i1 < NUMVARS; i1++) {
            first[j] = i0;
            second[j] = i1;
            combine[i0][i1] = j;
            combine[i1][i0] = j;
            j++;
        }
    }
    assert(j == NUMFCNS);
    
	/********************** done setting counters **********************/
    /*********************** DONE INITIALIZING ************************/
    
    /*********************************************************************************/
    /***********************START ITERATIVE OPTIMIZATION *****************************/
    /*********************************************************************************/
	
    Z_algor_frac_accum_algor_immed = 1;
    /* set lambdas to specified initial values */
	for(j=0; j< NUMFCNS; j++) {
		lambda[j] = init_lambda[j];
	}
   
	    
    /* initialize the sampler */
    initSampling();
    
    /* do Monte Carlo to get expectations with independent parameters */
    initSample(BURNINLENGTH, SEED);
	
    for(mc = 0; mc < NUMDATA; mc++) {
        sample();
    }
    
    /* done Monte Carlo to get expectations */
    
    /* start the clock */
    start_time = (double) clock(); 
    
    /* start iterating over new Monte Carlo runs */
    numruns = 0;
	/* add in exit condition */
	
    while((stop_me==0) && (numruns < NUMRUNS)) {
        numruns++;
 
        /* start (re)sampling */
        sample_time = (double) clock(); 
        initSample(BURNINLENGTH, SEED);
        
        for(j = 0; j < NUMFCNS; j++) {
            count_monte[j] = 0;
            bound_monte[j] = 1;
            sum_fjexpg_monte_eta[j] = 0;
            eta[j] = 0;
        }
        
        Z_monte_lambda = 0;
        Z_monte_eta = 0;
 
        for(mc = 0; mc < NUMSAMPS; mc++) {
            g_monte_lambda[mc] = sample();
            Z_monte_lambda += exp(g_monte_lambda[mc]);
            g_monte_eta[mc] = 0;
            Z_monte_eta++;
            
            if(sampleCountFcns > 0) {
                for(j = 0; j < sampleCountFcns; j++) {
                    i0 = prevSample[j];
                    if(count_monte[i0] >= bound_monte[i0]) {
                        bound_monte[i0] *= 2;
                        sInds_monte[i0] = (int *) realloc((void *) (sInds_monte[i0]), bound_monte[i0] * sizeof(int));
                        if(sInds_monte[i0] == NULL) {
                            fprintf(stderr, "Couldn't allocate to sInds_monte[%d], bound=%d, count=%d\n", i0, bound_monte[i0],count_monte[i0]);
                            exit(-1);
                        }
                    }
                    sInds_monte[i0][count_monte[i0]] = mc;
                    count_monte[i0]++;
                    sum_fjexpg_monte_eta[i0]++;
                }
            }
        }
    	for(j = 0; j < NUMFCNS; j++) {  
        	expec_algor[j] = ((double) count_monte[j] / (double) NUMSAMPS);
    	}
        sample_time = (double) clock()-sample_time; 
        
        /* start initializing monte checks */
        sum_gexpg_monte_eta = 0;
        sum_energyexpg_monte_eta = 0;
        /* done initializing monte checks */
        
        /* start learning from this Monte Carlo run */
        tsample = 0;
        while(tsample < TCOUNT) {
             
            
            /* start LEARNING ALGORITHM */
            learn_time = (double) clock();
            learn_alg();
            learn_time = (double) clock()-learn_time; 
        
            /* done LEARNING ALGORITHM */
            stop_me=WriteMetrics(trackFile,expec_empir,expec_algor,NUMVARS,numruns,tsample,start_time,sample_time, learn_time,mean_energy_monte_eta,rel_entropy_monte_theta_lambda,exit_conditions);
			sample_time=0;
            learn_time=0;
            
            tsample++;
        }
        /* done learning from this Monte Carlo run */
        
        /* start updating the lambdas for next run */
        for(j = 0; j < NUMFCNS; j++) {
            lambda[j] += eta[j];
        }
       
        /* done updating the lambdas for next run */
        
        
        
        Z_algor_frac_accum_algor_immed *= Z_monte_eta / NUMSAMPS;
		sprintf(numrunsString,"%d",numruns);
		strcpy(paramsFileTemp,paramsFile);
		strcat(paramsFileTemp,numrunsString);
		WriteParams(paramsFileTemp,lambda,first,second,NUMVARS,NUMFCNS);
        
        strcpy(paramsFileTemp,trackFile);
        strcat(paramsFileTemp,"Moments");
        strcat(paramsFileTemp,numrunsString);
        WriteParams(paramsFileTemp,expec_algor,first,second,NUMVARS,NUMFCNS);
        
    }
    /* done iterating over new Monte Carlo runs */
    /* DONE WITH ITERATIVE OPTIMIZATION */

    
	fflush(stderr);
    
    /* start comparing expectations */
    /* do Monte Carlo to get expectations with algorithm parameters */
	
    initSample(BURNINLENGTH, SEED);
    
    for(j = 0; j < NUMFCNS; j++) {
        expec_algor[j] = 0;
    }
	
   for(mc = 0; mc < NUMDATA; mc++) {
       for(j=0; j < NUMVARS; j++) {
           final_mc[j][mc]=0;
       }
   }
    
    for(mc = 0; mc < NUMDATA; mc++) {
        sample();
        
	
		if(sampleCountFcns > 0) {
            for(j = 0; j < sampleCountFcns; j++) {
                expec_algor[prevSample[j]]++;
					if(prevSample[j]<NUMVARS){
                    final_mc[prevSample[j]][mc] = 1;
                }
			}
        }
    }
    
    for(j = 0; j < NUMFCNS; j++) {
        expec_algor[j] /= NUMDATA;
	}
    fflush(stderr);

    /* end Monte Carlo to get expectations */
    /* start parameter output */
	WriteParams(paramsFile,lambda,first,second,NUMVARS,NUMFCNS);
	WriteMC(mcFile,final_mc,NUMVARS,NUMDATA);
	stop_me=WriteMetrics(trackFile,expec_empir,expec_algor,NUMVARS,numruns,tsample,start_time,sample_time, learn_time,mean_energy_monte_eta,rel_entropy_monte_theta_lambda,exit_conditions);
    /* DONE RESULTS */
    /* START FREEING MEMORY */
    fflush(stderr);
    
    endSampling();
	
    free(count_monte);
	/* free(init_lambda); */
    free(bound_monte);
    for(j = 0; j < NUMFCNS; j++) {
        free(sInds_monte[j]);
    }
    free(sInds_monte);
    free(lambda);
    free(g_monte_lambda);
    free(eta);
    free(g_monte_eta);
    free(sum_fjexpg_monte_eta);
    
    for(i0 = 0; i0 < NUMVARS; i0++) {
        free(sigma[i0]);
    }
    free(sigma);
    free(lambda_backup);
    free(expec_algor);
    free(first);
    free(second);
    for(i0 = 0; i0 < NUMVARS; i0++) {
        free(combine[i0]);
    }
    free(combine);
    free(final_mc);
    /* DONE FREEING MEMORY */
    
    return 0;
}
