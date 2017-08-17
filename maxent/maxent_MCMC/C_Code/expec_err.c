/**
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "learn_stop.h"
#include "expec_err.h"

#define INFO

void data_to_expec(int numvars, int numdata_s, int numdata_f,
                    int **data, double *expec) {
    int x, j, i1, i2;
    int numdata;
    double expecTot;
	numdata = numdata_f - numdata_s;
    
    expecTot = 0;
    for(j = 0; j < numvars; j++) {
      expec[j] = 0;
      for(x = numdata_s; x < numdata_f; x++) {
         expec[j] += data[j][x];
      }
      expec[j] /= numdata;
      expecTot += expec[j];
      /*fprintf(stderr, "-------expec[%d]: %f\n", j, expec[j]);*/
   }
   expecTot /= numvars;
   /*fprintf(stderr, "-------expecTot: %f\n", expecTot);*/
   
   for(i1 = 0; i1 < numvars; i1++) {
      for(i2 = i1 + 1; i2 < numvars; i2++) {
         expec[j] = 0;
         for(x = numdata_s; x < numdata_f; x++) {
            expec[j] += data[i1][x] * data[i2][x];
         }
         expec[j] /= numdata;
         j++;
      }
   }
}

void expec_diffs(int numvars, double *expecA, double *expecB,
                    double *avg_absdiff_expec1, double *avg_absdiff_expec2,
                    double *avg_absdiff_cov) {
    int j;
    int numfcns;
    
    double aad_e1, aad_e2, aad_cov;
    double connA, connB, var1, var2, covA, covB;
    double varridge; /*add a ridge to variances for numerical stability */
    
    
    #ifdef INFO
        double avg1_1, avg1_2, avg2_1, avg2_2;
    #endif
    
    numfcns = numvars + numvars * (numvars - 1) / 2;
	varridge= 0.0000000000000001;
		
    aad_e1 = 0;
    #ifdef INFO
        avg1_1 = 0; avg1_2 = 0; avg2_1 = 0; avg2_2 = 0;
    #endif
    for(j = 0; j < numvars; j++) {
        aad_e1 += fabs(expecA[j] - expecB[j]);
        #ifdef INFO
            avg1_1 += expecA[j]; avg1_2 += expecB[j];
        #endif
    }
    #ifdef INFO
        avg1_1 /= numvars; avg1_2 /= numvars;
    #endif
    aad_e2 = 0;
    aad_cov = 0;
    for(j = numvars; j < numfcns; j++) {
        connA = expecA[j] - expecA[first[j]] * expecA[second[j]];
        connB = expecB[j] - expecB[first[j]] * expecB[second[j]];
        
        var1 = expecA[first[j]] - expecA[first[j]] * expecA[first[j]];
        var2 = expecA[second[j]] - expecA[second[j]] * expecA[second[j]];
        covA = connA / (sqrt(var1)+varridge) / (sqrt(var2)+varridge);
        var1 = expecB[first[j]] - expecB[first[j]] * expecB[first[j]];
        var2 = expecB[second[j]] - expecB[second[j]] * expecB[second[j]];
        covB = connB / (sqrt(var1)+varridge) / (sqrt(var2)+varridge);
        aad_e2 += fabs(connA - connB);
        aad_cov += fabs(covA - covB);
        
        #ifdef INFO
            avg2_1 += expecA[j]; avg2_2 += expecB[j];
        #endif
    }
    #ifdef INFO
        avg2_1 /= (numfcns - numvars); avg2_2 /= (numfcns - numvars);
   /*     fprintf(stderr, "INFO: 1: %.14f, %.14f, 2: %.14f, %.14f\n", avg1_1, avg1_2, avg2_1, avg2_2);*/
    #endif
    
    aad_e1 /= numvars;
    aad_e2 /= (numfcns - numvars);
    aad_cov /= (numfcns - numvars);
    
    *avg_absdiff_expec1 = aad_e1;
    *avg_absdiff_expec2 = aad_e2;
    *avg_absdiff_cov = aad_cov;
}

/*int main(int argc, char *argv[]) {
    int j, x;
    int NUMVARS = 2;
    int NUMDATA = 5;
    int **sigma;
    double *output;
   
    sigma = (int **) malloc(NUMVARS * sizeof(int *));
    for(j = 0; j < NUMVARS; j++) {
        sigma[j] = (int *) malloc(NUMDATA * sizeof(int));
        for(x = 0; x < NUMDATA; x++) {
            sigma[j][x] = j + x;
        }
    }
    output = (double *) malloc(NUMVARS * sizeof(double *));
   
    data_to_expec(sigma, NUMVARS, NUMDATA, output);
    for(j = 0; j < NUMVARS; j++) {
        fprintf(stdout, "%f\n", output[j]);
    }

    for(j = 0; j < NUMVARS; j++) {
        free(sigma[j]);
    }
    free(sigma);

    return 0;
}   */
