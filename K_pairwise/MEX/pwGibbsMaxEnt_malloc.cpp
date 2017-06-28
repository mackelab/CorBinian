/*=================================================================
 * This is a MEX-file for MATLAB.  
 * All rights reserved.
 *=================================================================*/

/* pairwise Gibbs sampler MEX file */

#include <algorithm>   // for random_shuffle
#include "time.h"
#include "math.h"      // for exp()
namespace mat {
 #include "mex.h"   
}
//#include <vector>	   // for std::vector
//#include "random"    // for std::default_random_enginge


void mexFunction(int nlhs,mat::mxArray *plhs[],int nrhs,const mat::mxArray *prhs[]) {
    
/* 1. variable declarations */
    
    const int nSamples    = mat::mxGetScalar(prhs[0]); // number of samples 
    const int burnIn      = mat::mxGetScalar(prhs[1]); // length burnin phase
    const int d           = mat::mxGetScalar(prhs[2]); // data dimensionality
    	  double * x0          = mat::mxGetPr(prhs[3]); // initial chain segment
    	  double * pairs0      =  mat::mxGetPr(prhs[4]); // look-up table #1
    	  double * fm0         =  mat::mxGetPr(prhs[6]); // look-up table #2
    const double * h           = mat::mxGetPr(prhs[7]); // single-cell weights
    const double * J           = mat::mxGetPr(prhs[8]); // cell-pair weights
    const double * L           = mat::mxGetPr(prhs[9]); // population weights
    
    int numAll  = d*(d+3)/2+1;	// number of features for (h,J,V(K)) model
    int numPair = d*(d-1)/2;    // number of pairs between d variables 
	int numIsing = numPair+d;   // number of features for Ising model 
   	 
    int * pairs            = new int[2*numPair];
    int * fm               = new int[d*(d-1)];
    
    // The following variables are of dynamic size (=> quite the pain in C)
    mat::mxArray *outVar = mat::mxCreateDoubleMatrix(numAll, 1, mat::mxREAL);
    double * xSampled = mat::mxGetPr(outVar);
//     double * xTmp = malloc(numAll * sizeof(double)); // intermediate storage (results after each sweep)
//     double * pc   = malloc(numAll * sizeof(double)); // short-term storage (results within each sweep)
//     bool   * xc   = malloc(d * sizeof(bool)); // current chain member   
//     int    * ks   = malloc(numPair * sizeof(int)); // index vector to keep track of pairs in a sweep     
    double * xTmp = new double[numAll]; // intermediate storage (results after each sweep)
    double * pc   = new double[numAll]; // short-term storage (results within each sweep)
    bool   * xc   = new bool[numAll]; // current chain member   
    
    int    * ks   = new int[numPair]; // index vector to keep track of pairs in a sweep 
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
    
    srand((unsigned)time(NULL));              // seed random stream

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
    
    for (idx=0; idx<2*numPair; idx++) {
        *(pairs+idx) = (int) *(pairs0+idx);  // 
    }       
    
    for (idx=0; idx<d*(d-1); idx++) {
        *(fm+idx) = (int) *(fm0+idx);  // 
    }       
    
    for (idx=0; idx<numAll; idx++) {
        xTmp[idx] = 0;  // intermediate storage for Rao-Blackwell variables
    }       
/* 3. Start MCMC chain */
 
	/* FOR EACH OF THE IN TOTAL nSamples + burnIn MANY SWEEPS... */
    for (i = -burnIn+1; i < nSamples+1; i++) {

    	std::random_shuffle(ksfirst, kslast); // shuffle order for current sweep
        
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
                            
        } // end for j = 0:numPair

        /* normalize results by number of visits */
        for (idx=0; idx<d; idx++) {
        	pc[idx] /= (d-1);
            pc[numIsing+idx] /= numPair;
        }
        pc[numAll-1] /= numPair; // final entry for L (d+1 entries!)
        
        /* shove results to mid-term storage and final outputs */
		if (i > 0) {          //* after burn-in phase is over:
        	for (idx =0; idx<numAll; idx++) {// copy into mid-term storage
        		xTmp[idx] += pc[idx]; // adds up results from 100 sweeps 
            }                         // before passing on to xSampled
            
        	if ((i % 100) == 0) {    // every thousand samples ...
            	for (idx = 0; idx<numAll; idx++) { // move to final results 
            		xSampled[idx] = ((i-100)*xSampled[idx] + xTmp[idx])/i;	
            		xTmp[idx] = 0; // (normalize regularly to avoid 
            	}                  //  numerical overflows)
                //mexPrintf("\nCurrent sample is # %d ", i);   
                //mexPrintf(" out of %d .", nSamples);   
        	}
        }
        else { // give status reports also during burn-in phase
            if (((i+burnIn) % 100) == 0) {
                //mexPrintf("\nAt burnin step # %d ", (i+burnIn));   
                //mexPrintf(" out of %d .", burnIn);               
            }
        }
    } // end for i = -burnIn:nSamples 
    
 plhs[0] = outVar;
 if (nlhs > 1) {
    mat::mxArray *outX = mat::mxCreateDoubleMatrix(d, 1, mat::mxREAL);
    double * xOut = mat::mxGetPr(outX);
    for (idx=0; idx<d; idx++) {
        xOut[idx]=0;
        if (xc[idx]) {
            xOut[idx] = 1;
        }
    }
    plhs[1] = outX;    
 }
 delete [] pairs;
 delete [] fm;
 delete [] xTmp;
 delete [] pc;
 delete [] xc;
 delete [] ks;
 
 //plhs[2] = outE;  
}    

