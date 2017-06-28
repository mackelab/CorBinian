/**
 * Implements a Gibbs sampler.
 *
 * COMPILE: make sampler
 *
 * RUN: sampler.x
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "mt19937ar.h"
#include "sampler.h"
#include "learn_stop.h"

void initSampling() {
    prevSample = (int *) malloc(NUMFCNS * sizeof(int));
    assert(prevSample != NULL);
}

void initSample(int burnin_length, unsigned long int seed) {
    int i;
    
    assert(prevSample != NULL);
    
    /* initialize the random number generator */
    init_genrand(seed);
    
    /* initialize the sample */
    sampleCountVars = 0;
    sampleCountFcns = 0;
    
    for(i = 0; i < burnin_length; i++) {
        sample();
    }
}

void endSampling() {
    free(prevSample);
}

/* Returns g_{\lambda} under the current distribution.
 * g_{\lambda_temp} under the new distribution, which will be constructed during the
 * following round, is always 0 since each \lambda_temp is initialized to 0
 */
double sample() {
    int n, i, j, k, isOne, eqVar;
    double p, r;
    
    for(n = 0; n < NUMVARS; n++) {
        i = floor(NUMVARS * genrand_real2());
        p = 0;
        isOne = 0;
        for(j = 0; j < sampleCountVars; j++) {
            k = prevSample[j];
            if(k == i) {
                isOne = 1;
                eqVar = j;
            } else {
                p += lambda[combine[i][k]];
            }
        }
        p = 1.0/(1.0 + exp(-lambda[i] - p));

        r = genrand_real2();
        if(r < p) {
            /* set the i_th coordinate to 1 if it is not already 1 */
            if(!isOne) {
                prevSample[sampleCountVars] = i;
                sampleCountVars++;
            }
        } else if(isOne) {
            for(j = eqVar; j < sampleCountVars - 1; j++) { prevSample[j] = prevSample[j+1]; }
            sampleCountVars--;
        }
    }
    
    return initLI();
}

double initLI() {
    double Dbl;
    
    int i, j, k, l;

    Dbl = 0;
    sampleCountFcns = sampleCountVars;
    if(sampleCountVars == 0) return Dbl;
        
    for(i = 0; i < sampleCountVars; i++) {
        k = prevSample[i];
        Dbl += lambda[k];
        
        for(j = i+1; j < sampleCountVars; j++) {
            l = prevSample[j];
            prevSample[sampleCountFcns] = combine[k][l];
            sampleCountFcns++;
            Dbl += lambda[combine[k][l]];
        }
    }
    
    return Dbl;
}
        
/* a main for testing */
/* needs to be updated before use */
/*int main(int argc, char *argv[]) {
    int n, i;
    int counts[10];
    
    for(i = 0; i < 10; i++) {
        counts[i] = 0;
    }
    
    init_genrand(4);
    for(n = 0; n < 10000000; n++) {
        i = floor(genrand_real2() * 10);
        counts[i]++;
    }
    
    for(i = 0; i < 10; i++) {
        fprintf(stdout, "%d ", counts[i]);
    }
    fprintf(stdout, "\n");
    
    return 0;
}*/

