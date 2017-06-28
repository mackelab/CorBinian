/************************************************************************
 * Schapire's maximum entropy algorithm
 ************************************************************************/
 
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "learn_iterative.h"
#include "learn_stop.h"
#define max(A,B) ((A)>(B)) ? (A) : (B)
#define min(A,B) ((A)<(B)) ? (A) : (B)


void learn_alg() {
    int minJ;
    double minDel;
    
    int i0, j, mc;
    double estimfj, expdelta, DelL, mindelta;
    double oldg, newg; 
    double minallowed, minstep;
    /* START OPTIMIZING PARAMETERS */
    /* find the best lambda to change */
    minDel = 0;
    minJ = 0;
	minallowed=0.0000001;
	minstep=exp(-20);
	
    for(j = 0; j < NUMFCNS; j++) {
        estimfj = sum_fjexpg_monte_eta[j] / Z_monte_eta;
		/* make sure that we do not get 0 or 1s in the ratio, for numerical stability */
		/*of course, it would be way more efficient to do these checks only once, not in every inner loop... */
		/*estimfj=max(estimfj,minallowed);
		estimfj=min(estimfj,1-minallowed);	
		expec_empir[j]=max(expec_empir[j],minallowed);
		expec_empir[j]=min(expec_empir[j],1-minallowed);*/
		
        expdelta = expec_empir[j] * (1 - estimfj);
        expdelta /= (1 - expec_empir[j]) * estimfj;
		
		/*expdelta=min(expdelta,maxstep);
		expdelta=max(expdelta,minstep); */
		
		DelL = -log(expdelta) * expec_empir[j];
        DelL += log(1 + (expdelta - 1) * estimfj);
        
		
        if(DelL < minDel && (expdelta>minstep) && (1 + (expdelta - 1) * estimfj>minstep)) {
            minDel = DelL;
            minJ = j;
            mindelta = log(expdelta);
        }
    }
	/*for numerical stability, limit maximal stepsize we allow: */
	mindelta=max(mindelta,-5);
	mindelta=min(mindelta,5);
	
    /*fprintf(stderr, "exp[minJ]: %f, mindelta: %f, minJ: %d\n", eta[minJ], mindelta, minJ);
    */
    /* DONE OPTIMIZING PARAMETERS */
    
    /* START UPDATING */
    eta[minJ] += mindelta;
    /* start Monte Carlo run updates */
    /*fprintf(stderr, "Z0: %f\n", Z_monte_eta);*/
    for(i0 = 0; i0 < count_monte[minJ]; i0++) {
        mc = sInds_monte[minJ][i0];
        oldg = g_monte_eta[mc];
        g_monte_eta[mc] += mindelta;
        newg = g_monte_eta[mc];
        
        Z_monte_eta += exp(newg) - exp(oldg);
        
        /* start monte checks */
        sum_gexpg_monte_eta += newg * exp(newg) - oldg * exp(oldg);
        sum_energyexpg_monte_eta += g_monte_lambda[mc] * (exp(newg) - exp(oldg));
        /* done monte checks */
    } 
    /*fprintf(stderr, "Z1: %f\n", Z_monte_eta);*/
    /* done Monte Carlo run updates */
    
    /* start function expectation updates */
    for(j = 0; j < NUMFCNS; j++) {
        sum_fjexpg_monte_eta[j] = 0;
        for(i0 = 0; i0 < count_monte[j]; i0++) {
            sum_fjexpg_monte_eta[j] += exp(g_monte_eta[sInds_monte[j][i0]]);
        }
    }
    /* done function expectation udpates */
    
    /* start final checks */
    /* start monte */
    rel_entropy_monte_theta_lambda = -log(Z_monte_eta / NUMSAMPS) + sum_gexpg_monte_eta / Z_monte_eta;
    mean_energy_monte_eta = sum_energyexpg_monte_eta / Z_monte_eta;
    /* done monte */
    
    /* done final checks */
    /* DONE UPDATING */
        
    return;
}
