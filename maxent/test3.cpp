/*=================================================================
 * mexfunction.c  
 * This is a MEX-file for MATLAB.  
 * All rights reserved.
 *=================================================================*/

#include "mex.h"
#include "time.h"
//#include <iostream>  // for cout
#include "math.h"      // for exp()
//#include <algorithm> // for random_shuffle
//#include <vector>	   // for std::vector
//#include "random"    // for std::default_random_enginge


void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    /* 1. variable declarations */
    const int nSamples    = mxGetScalar(prhs[0]); // number of samples 
    const int burnIn      = mxGetScalar(prhs[1]); // length burnin phase
    const int d           = mxGetScalar(prhs[2]); // data dimensionality
    const double * x0          = mxGetPr(prhs[3]); // initial chain segment
    const double * pairs       = mxGetPr(prhs[4]); // look-up table #1
    const double * fm          = mxGetPr(prhs[6]); // look-up table #2
    const double * h           = mxGetPr(prhs[7]); // single-cell weights
    const double * J           = mxGetPr(prhs[8]); // cell-pair weights
    const double * L           = mxGetPr(prhs[9]); // population weights
    
    int numAll  = d*(d+3)/2+1;	// number of features for (h,J,V(K)) model
    int numPair = d*(d-1)/2;    // number of pairs between d variables 
	int numIsing = numPair+d;   // number of features for Ising model 
   	
	// The following variables are of dynamic size (=> quite the pain in C)
    mxArray *outVar = mxCreateDoubleMatrix(136, 1, mxREAL);
    double * xSampled = mxGetPr(outVar);
    double   xTmp[136]; // intermediate storage (results after each sweep)
    double     pc[136]; // short-term storage (results within each sweep)
    bool        xc[15]; // current chain member   
    
    int        ks[105]; // index vector to keep track of pairs in a sweep 
    int* ksfirst = ks;           // pointers to first and last element of
    int* kslast  = (ks+numPair); // index array
    
    int idl=0;					// current population activity count K
    
    int k,l;					// indices for current pair of variables 
    int dmok, dmol, fmidx;      // convenience indexes for entries of J
    double sk, sl, slk;			// summary variable for effects from J
    double L0, L1, L2;		    // variables for effects from L

    double p[4]; // probabilities for bivariate distribution over x(k),x(l)
    double sump, maxp; // normalizer and 
  
    double rnd;	// random variable to decide new entries for x(k), x(l)

    int idx; // all-purpose index for iterating through vector entries
    int i,j; // indexes for outer and inner Gibbs loops, respectively
    
    /* 2. Initialize variables where needed */
//  srand((unsigned)time(NULL));              // seed random stream

    for (idx=0; idx<numPair; idx++) {
        ks[idx] = idx; // default order of pairs: counted from 1 to numPair
    }
    for (idx=0; idx<d; idx++) {
    	if (x0[idx]==1) {
            xc[idx] = true; // current chain segment x
    		idl+=1;	        // population activity count
        }
        else {
            xc[idx] = false;
    	}
    }
    
    for (idx=0; idx<numAll; idx++) {
        xTmp[idx] = 0;  // intermediate storage for Rao-Blackwell variables
    }
        
	/* FOR EACH OF THE IN TOTAL nSamples + burnIn MANY SWEEPS... */
    for (i = -burnIn+1; i < nSamples+1; i++) {

//     	random_shuffle(ksfirst, kslast); // shuffle order for current sweep
        
        for (idx=0; idx < numAll; idx++) {
        	pc[idx] = 0;	// reset current results for expected values
        }

        /* ... VISIT EACH POSSIBLE PAIR OF VARIABLES AND UPDATE THEM */
        for (j = 0; j < numPair; j++) {
     
            k = *(pairs+ks[j]);			// choose current pair of 
         	l = *(pairs+ks[j]+numPair);	// variables to work on
            if (xc[k]) { //
                idl--;   // momentarily assume both variables
            }            // to equal 0. 
            if (xc[l]) { // We'll add the true values at the 
                idl--;   // end of the loop iteration.
            }            //
            
            /* compute probabilities for bivariate distribution over x(k),
             * x(l) by inspecting all four possible configurations */
    		xc[k] = true;
    		xc[l] = false;
    		sk  = 0; // storage for effects from J for case x(k)=1, x(l)=0
            dmok  = (d-1)*k; // column index in vectorized look-up table fm
    		for (idx=0; idx < d-1; idx++) {
                fmidx = *(fm+idx+dmok); // current index within pairs 
                if (xc[(int) *(pairs+fmidx)] && xc[(int) *(pairs+fmidx+numPair)]) {
                    sk += J[fmidx];
                }
            }
    		xc[k] = false;
    		xc[l] = true;
    		sl  = 0; // storage for effects from J for case x(k)=0, x(l)=1
            dmol  = (d-1)*l; // column index in vectorized look-up table fm
    		for (idx=0; idx < d-1; idx++) {
                fmidx = *(fm+idx+dmol); // current index within pairs 
                if (xc[(int) *(pairs+fmidx)] && xc[(int) *(pairs+fmidx+numPair)]) {
                    sl += J[fmidx];
                    
                }
            }
    		xc[k] = true; // if both x(k),x(l)=1, we need to add a few more
    		slk = 0;      // entries of J that will be stored in slk
    		for (idx=0; idx < d-1; idx++) {
                fmidx = *(fm+idx+dmok); // current index within pairs 
                if (xc[(int) *(pairs+fmidx)] && xc[(int) *(pairs+fmidx+numPair)]) {
                    slk += J[fmidx];
                }
            }
            L0 = *(L+idl);   // retrieve entries of L corresponding
    		L1 = *(L+idl+1); // to 0,1 or 2 out of the variables
    		L2 = *(L+idl+2); // x(k) and x(l) being equal to 1

    		/* finish computing bivariate Bernoulli probabilities */
    		p[0]= (h[k]        + sk        + L1); // x(k) = 1, x(l) = 0
    		p[1]= (h[k] + h[l] + slk + sl  + L2); // x(k) = 1, x(l) = 1
    		p[2]= (       h[l]       + sl  + L1); // x(k) = 0, x(l) = 1
    		p[3]= (                          L0); // x(k) = 0, x(l) = 0
            
            maxp = 0; // maximum of p[0],...,p[3]
            for (idx=0; idx<4; idx++) {
             if (p[idx]>maxp) {
                 maxp = p[idx];
             }
            }
            for (idx=0; idx<4; idx++) { // subtract maxp to hopefully 
             p[idx] = exp(p[idx]-maxp); // avoid numerical overflows
            }
    		sump = p[0]+p[1]+p[2]+p[3];
            p[0] /= sump;									//
            p[1] /= sump;									// normalize
            p[2] /= sump;									//
            p[3] /= sump;									//
            
            /* draw new pair of variables from bivariate Bernoulli */
            rnd = (rand()+1.0)/(RAND_MAX+1.0); // c99-style uniform distr.
            if (rnd>(p[0]+p[1])) { //
                xc[k] = false;     // we last visited the configuration 
            }                      // x(k)=x(l)=1, now just need to flip
            else {                 // x(k),x(l) to zero if necessary. 
                idl++;             //
            }
            if ((p[0] > rnd) || (rnd > (1 - p[3]))) {
                xc[l] = false;
            }
            else {
                idl++;
            }

            /* collect probabilities for Rao-Blackwell */
            pc[k] += (p[0]+p[1]);
            pc[l] += (p[1]+p[2]);
            pc[d+ks[j]] = p[1];
        	pc[numIsing+idl]++;
            
            
//             if (xc[0]) {
//                 mexPrintf("\n x(0) = ( %d ", 1);
//             }
//             else {
//                 mexPrintf("\n x(0) = ( %d ", 0);
//             }
//                 mexPrintf("), ");
//             if (xc[1]) {
//                 mexPrintf(", x(1) = ( %d ", 1);
//             }
//             else {
//                 mexPrintf(", x(1) = ( %d ", 0);
//             }                mexPrintf("), ");
//                 mexPrintf("), ");
//             if (xc[2]) {
//                 mexPrintf(", x(2) = ( %d ", 1);
//             }
//             else {
//                 mexPrintf(", x(2) = ( %d ", 0);
//             }
//                 mexPrintf("), ");
//             if (xc[3]) {
//                 mexPrintf(", x(3) = ( %d ", 1);
//             }
//             else {
//                 mexPrintf(", x(3) = ( %d ", 0);
//             }                
//                 mexPrintf("), ");
//                 mexPrintf(", idl = %d ", idl);
                
        } // end for j = 0:numPair

        /* normalize results by number of visits */
        for (idx=0; idx<d; idx++) {
        	pc[idx] /= (d-1);
            pc[numIsing+idx] /= numPair;
        }
        pc[numAll-1] /= numPair; // final entry for L; d+1 entries for population activity count since 0 <= idl <= d
        for (idx = 0; idx < numAll; idx++) {
//             mexPrintf("\npc( %d ", idx+1);  
//             mexPrintf(") is %f .", pc[idx]);  
        } 
//         mexPrintf("\n----------------------------------");  
        
		if (i > 0) {          //* after burn-in phase is over:
        	for (idx = 0; idx < numAll; idx++) {					// copy into mid-term storage
        		xTmp[idx] += pc[idx]; // adds up results from 1000 sweeps before passing on to xSampled
//                     mexPrintf("\nxTmp( %d ", idx+1);  
//                     mexPrintf(") is %f .", xTmp[idx]);  
            }
//             mexPrintf("\n*************************************");  
            
        	if ((i % 1000) == 0) {										// every thousand samples ...
            	for (idx = 0; idx < numAll; idx++) {
            		xSampled[idx] = ((i-1000)*xSampled[idx] + xTmp[idx])/i;	// ... move to final results (normalize
            		xTmp[idx] = 0;                                        // regularly to avoid numerical overflows)
            	}
                mexPrintf("\nCurrent sample is %d ", i);   
                mexPrintf(" out of %d .", nSamples);   
        	}
        }
        else {
            if (((i+burnIn) % 100) == 0) {
                mexPrintf("\nAt burnin step is %d ", (i+burnIn));               
            }
        }
        
    } // end for i = -burnIn:nSamples 

//      for (idx = 0; idx < numAll; idx++) {
//         mexPrintf("\nxSampled( %d ", idx+1);  
//         mexPrintf(") is %f .", xSampled[idx]);   
//      } 
//      for (idx = 0; idx < numPair; idx++) {
//         tmp = *(pairs+idx);  
//         mexPrintf("\n pairs(:,1) %f .", tmp);   
//      } 
//      for (idx = 0; idx < numPair; idx++) {
//         tmp = *(pairs+idx+numPair);  
//         mexPrintf("\n pairs(:,2) %f .", tmp);   
//      } 
//             mexPrintf("\n fm:");
//             for (idx=0; idx<(d*(d-1)); idx++) {
//                 tmp = *(fm+idx);
//                 mexPrintf("\n %d", tmp);
//             }
//             mexPrintf("\n pairs:");
//             for (idx=0; idx<2*numPair; idx++) {
//                 tmp = *(pairs+idx);
//                 mexPrintf("\n %d", tmp);
//             }
    
   plhs[0] = outVar;
}    

