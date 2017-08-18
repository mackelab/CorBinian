/**
 * A header file for sampler algorithms.
 */
 
#ifndef SAMPLER_H
#define SAMPLER_H

int *prevSample;

int sampleCountVars;
int sampleCountFcns;

void initSampling();
void endSampling();
void initSample(int burnin_length, unsigned long int seed);
double sample();
double initLI();

#endif
