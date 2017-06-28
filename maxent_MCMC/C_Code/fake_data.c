/**
 * Run a Monte Carlo to generate "fake" data from an
 * input file of h's and J's
 *
 * COMPILE: make fake
 *
 * RUN: fake_data.x -v 10 -d 200000 -b 200000 -D ../../DATA/n40.txt > ../../DATA/fake_data.v10.d200000.b200000.out
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "learn_stop.h"
#include "sampler.h"

#define SEED 1

#define TRUE 1
#define FALSE 0

int main(int argc, char *argv[]) {
    /* START DEFINING */
    int BURNINLENGTH;
    
    /* start counters, miscellaneous */
    int a, j, i0, i1;
    double dtmp;
    int REQ;
    FILE *fhj;
    char dataFile[200];
    /* done counters, miscellaneous */
    
    /* start checking */
    double Z, H;
    double *expec;
    /* done checking */
    /* DONE DEFINING */
    
    /* START INITIALIZING */
    /* start reading inputs */
    REQ = 0;
    for(a = 1; a < argc; a++) {
        if(argv[a][0] == '-') {
            switch(argv[a][1]) {
                case 'v':
                    NUMVARS = atoi(argv[a] + 3); REQ ^= (1<<0); break;
                case 'd':
                    NUMDATA = atoi(argv[a] + 3); REQ ^= (1<<1); break;
                case 'b':
                    BURNINLENGTH = atoi(argv[a] + 3); REQ ^= (1<<2); break;    
                case 'D':
                    strcpy(dataFile, argv[a] + 3); REQ ^= (1<<3); break;
                default:
                    fprintf(stderr, "Bad option: %s\n", argv[a]); exit(-1);
            }
        }
    }
    if(REQ != (1<<4) - 1) {
        fprintf(stderr, "Usage: %s\n", argv[0]);
        /****** more info on Usage ******/
        exit(-1);
    }
    /* done reading inputs */
    
    /* start allocating space */
    NUMFCNS = NUMVARS + NUMVARS * (NUMVARS - 1) / 2;
    /* h's and J's */
    lambda = (double *) malloc(NUMFCNS * sizeof(double));
    /* start samples */
    sigma = (int **) malloc(NUMVARS * sizeof(int *));
    for(i0 = 0; i0 < NUMVARS; i0++) {
        sigma[i0] = (int *) malloc(NUMDATA * sizeof(int));
    }
    /* done samples */
    
    /* start counters */
    first = (int *) malloc(NUMFCNS * sizeof(int));
    second = (int *) malloc(NUMFCNS * sizeof(int));
    combine = (int **) malloc(NUMVARS * sizeof(int *));
    for(i0 = 0; i0 < NUMVARS; i0++) {
        combine[i0] = (int *) malloc(NUMVARS * sizeof(int));
    }
    /* done counters */
    
    /* start checking */
    expec = (double *) malloc(NUMFCNS * sizeof(double));
    /* done allocating space */
    
    /* start setting counters */
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
    /* done setting counters */
    /* DONE INITIALIZING */

    /* START READING DATA */
    if((fhj = fopen(dataFile, "r")) == NULL) {
        fprintf(stderr, "Bad file name: %s\n", dataFile);
        exit(-1);
    }
    
    j = 0;
    a = TRUE;
    do {
        if(fscanf(fhj, "%d %lf", &i0, &dtmp) == EOF) {
            fprintf(stderr, "Couldn't read h's properly.\n");
            exit(-1);
        }
        
        if(i0 == j) {
            lambda[j] = dtmp;
        } else {
            a = FALSE;
        }
        j++;
    } while((a == TRUE) && (j < NUMVARS));
    if(j != NUMVARS) {
        fprintf(stderr, "Too few h's.\n");
        exit(-1);
    }
    while(fscanf(fhj, "%d %d %lf", &i0, &i1, &dtmp) != EOF) {
        if((i0 == first[j]) && (i1 == second[j])) {
            lambda[j] = dtmp;
            j++;
        }
    }
    if(j != NUMFCNS) {
        fprintf(stderr, "Wrong number of J's.\n");
        exit(-1);
    }
    fclose(fhj);
    /* DONE READING DATA */
    
    /* START SAMPLING */
    /* initialize the sampler */
    initSampling();
    
    initSample(BURNINLENGTH, SEED);
    for(a = 0; a < NUMDATA; a++) {
        sample();
        
        /* fill in the data array for this sample */
        for(i0 = 0; i0 < NUMVARS; i0++) {
            sigma[i0][a] = 0;
        }
        for(i0 = 0; i0 < sampleCountVars; i0++) {
            sigma[prevSample[i0]][a] = 1;
        }
    }
    /* DONE SAMPLING */
    
    /* START PRINTING */
    for(i0 = 0; i0 < NUMVARS; i0++) {
        for(a = 0; a < NUMDATA - 1; a++) {
            fprintf(stdout, "%d ", sigma[i0][a]);
        }
        fprintf(stdout, "%d\n", sigma[i0][a]);
    }
    /* DONE PRINTING */
    
    /* START CHECKING */
    Z = 0;
    for(j = 0; j < NUMVARS; j++) {
        expec[j] = 0;
    }
    for(a = 0; a < (1 << NUMVARS); a++) {
        H = 0;
        for(i0 = 0; i0 < NUMVARS; i0++) {
            H += ((a >> i0)&1) * lambda[i0];
            for(i1 = i0 + 1; i1 < NUMVARS; i1++) {                           
                H += ((a >> i0)&1) * ((a >> i1)&1) * lambda[combine[i0][i1]];
            }
        }   
        
        for(j = 0; j < NUMVARS; j++) {
            expec[j] += ((a >> j)&1) * exp(H);
        }
        for( ; j < NUMFCNS; j++) {
            expec[j] += ((a >> first[j])&1) * ((a >> second[j])&1) * exp(H);
        }
        
        Z += exp(H);
    }
    
    for(j = 0; j < NUMFCNS; j++) {
        fprintf(stderr, "%d %f\n", j, expec[j]/Z);
    }
    /* DONE CHECKING */
    
    /* START FREEING MEMORY */
    endSampling();
    
    free(lambda);
    for(i0 = 0; i0 < NUMVARS; i0++) {
        free(sigma[i0]);
    }
    free(sigma);
    free(first);
    free(second);
    for(i0 = 0; i0 < NUMVARS; i0++) {
        free(combine[i0]);
    }
    free(combine);
    /* DONE FREEING MEMORY */

    return 0;
}
