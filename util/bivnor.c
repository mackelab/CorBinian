/*
 * =============================================================
 * binor.c
 *
 * Approximation of the cumulative bivariate normal probability
 *
 * This is a MEX-file for MATLAB.
 * Copyright (c) 2004 William Moranvil (moranviw@onid.orst.edu)
 * =============================================================
 */
 
/* $Revision: 1.9 $ */

#include "mex.h"
#include "math.h"

#define SQRT2 1.414213562373095049
#define SQRTPI 1.772453850905516027
#define twopi 6.283195307179587

/* ---------------------------------------------------------------------------
                    UNIVARIATE NORMAL PROBABILITY
   
   Compute the cumulative distribution function of N(0,1), ie normal with 
   zero mean and unit variance.
   The C code was written by Ajay Shah (ajayshah@usc.edu)
   I only had to change few stuff for the code to be Matlab C compatible
   ---------------------------------------------------------------------------*/
#define UPPERLIMIT 20.0 
/* will not return either of univariate normal density or 
	probability when x < -UPPERLIMIT  or x > UPPERLIMIT. */

#define P10 242.66795523053175
#define P11 21.979261618294152
#define P12 6.9963834886191355
#define P13 -.035609843701815385
#define Q10 215.05887586986120
#define Q11 91.164905404514901
#define Q12 15.082797630407787
#define Q13 1.0

#define P20 300.4592610201616005
#define P21 451.9189537118729422
#define P22 339.3208167343436870
#define P23 152.9892850469404039
#define P24 43.16222722205673530
#define P25 7.211758250883093659
#define P26 .5641955174789739711
#define P27 -.0000001368648573827167067
#define Q20 300.4592609569832933
#define Q21 790.9509253278980272
#define Q22 931.3540948506096211
#define Q23 638.9802644656311665
#define Q24 277.5854447439876434
#define Q25 77.00015293522947295
#define Q26 12.78272731962942351
#define Q27 1.0

#define P30 -.00299610707703542174
#define P31 -.0494730910623250734
#define P32 -.226956593539686930
#define P33 -.278661308609647788
#define P34 -.0223192459734184686
#define Q30 .0106209230528467918
#define Q31 .191308926107829841
#define Q32 1.05167510706793207
#define Q33 1.98733201817135256
#define Q34 1.0

double pnorm1(double phi[], double x[])
{
	int sn;
	double R1, R2, R3, y, y2, y3, y4, y5, y6, y7;
	double erf, erfc, z, z2, z3, z4;
	/*double phi;*/

	if (x[0] < -UPPERLIMIT){
	phi[0]=0.0;
	return;
	}
	
	if (x[0] > UPPERLIMIT){
	phi[0]=1.0;
	return;
	}

	y = x[0] / SQRT2;
	if (y < 0) {
		y = -y;
		sn = -1;
	}
	else
		sn = 1;

	y2 = y * y;
	y4 = y2 * y2;
	y6 = y4 * y2;

	if(y < 0.46875) {
		R1 = P10 + P11 * y2 + P12 * y4 + P13 * y6;
		R2 = Q10 + Q11 * y2 + Q12 * y4 + Q13 * y6;
		erf = y * R1 / R2;
		if (sn == 1)
			phi[0] = 0.5 + 0.5*erf;
		else 
			phi[0] = 0.5 - 0.5*erf;
	}
	else 
		if (y < 4.0) {
			y3 = y2 * y;
			y5 = y4 * y;
			y7 = y6 * y;
			R1 = P20 + P21 * y + P22 * y2 + P23 * y3 + 
			    P24 * y4 + P25 * y5 + P26 * y6 + P27 * y7;
			R2 = Q20 + Q21 * y + Q22 * y2 + Q23 * y3 + 
			    Q24 * y4 + Q25 * y5 + Q26 * y6 + Q27 * y7;
			erfc = exp(-y2) * R1 / R2;
			if (sn == 1)
				phi[0] = 1.0 - 0.5*erfc;
			else
				phi[0] = 0.5*erfc;
		}
		else {
			z = y4;
			z2 = z * z;
			z3 = z2 * z;
			z4 = z2 * z2;
			R1 = P30 + P31 * z + P32 * z2 + P33 * z3 + P34 * z4;
			R2 = Q30 + Q31 * z + Q32 * z2 + Q33 * z3 + Q34 * z4;
			erfc = (exp(-y2)/y) * (1.0 / SQRTPI + R1 / (R2 * y2));
			if (sn == 1)
				phi[0] = 1.0 - 0.5*erfc;
			else 
				phi[0] = 0.5*erfc;
		}
	return;
}

/*---------------------------------------------------------------------------

                    BIVARIATE NORMAL PROBABILITY
    
  based on algorithm 462 Commun. ACM oct 1973 p638
  gives the probability that a bivariate normal exceeds (ah,ak).
  gh and gk are .5 times the right tail areas of ah, ak under a n(0,1)

  Tranlated from FORTRAN to ratfor using struct; from ratfor to C by hand.
  The C code was written by Ajay Shah (ajayshah@usc.edu)
  I only had to change few stuff for the code to be Matlab C compatible

  ---------------------------------------------------------------------------*/
  
#define con (twopi / 2.0) * 10.0e-10

double bivnor(double output[], double ah[], double ak[], double r[])
{
	double a2, ap, b, cn, conex, ex, g2, gh, gk, gw, h2, h4, rr, s1, s2, 
	sgn, sn, sp, sqr, t, temp, w2, wh, wk;
	int is;
	
	/* necessary to call subroutine pnorm1*/
	double *buf1;					
	double *buf2;
	
    /* memory allocation */
    if ((buf1=(double*) mxMalloc(1 * sizeof(double)))==NULL)   
        mexErrMsgTxt("l.51 memory allocation failed.");
    if ((buf2=(double*) mxMalloc(1 * sizeof(double)))==NULL)   
        mexErrMsgTxt("l.51 memory allocation failed.");  

	buf1[0] = -ah[0];
	pnorm1(buf2,buf1);
	gh = buf2[0];
	gh = gh / 2.0;
	
	buf1[0] = -ak[0];
	pnorm1(buf2,buf1);
	gk = buf2[0];
	gk = gk / 2.0;

	b = 0;

	if (r[0]==0)
		b = 4*gh*gk;
	else {
		rr = 1-r[0]*r[0];
		
		if (rr<0){                                                                     /* Error if 1-r^2 < 0 */
		    output[0]=-99.9;
		    /* free all memory allocated dinamically before terminate routine*/
	        mxFree(buf1);
	        mxFree(buf2);
			return;}
			
		if (rr!=0) {
			sqr = sqrt(rr);
			if (ah[0]!=0) {
				b = gh;
				if (ah[0]*ak[0]<0)
					b = b-.5;
				else if (ah[0]*ak[0]==0)
					goto label10;
			}
			else if (ak[0]==0) {
				b = atan(r[0]/sqr)/twopi+.25;
				goto label50;
			}
			b = b+gk;
			if (ah[0]==0)
				goto label20;
label10:
			wh = -ah[0];
			wk = (ak[0]/ah[0]-r[0])/sqr;
			gw = 2*gh;
			is = -1;
			goto label30;
label20:
			do {
				wh = -ak[0];
				wk = (ah[0]/ak[0]-r[0])/sqr;
				gw = 2*gk;
				is = 1;
label30:
				sgn = -1;
				t = 0;
				if (wk!=0) {
					if (fabs(wk)>=1)
						if (fabs(wk)==1) {
							t = wk*gw*(1-gw)/2;
							goto label40;
						}
						else {
							sgn = -sgn;
							wh = wh*wk;
							buf1[0]=wh;
							pnorm1(buf2,buf1);
							g2 =buf2[0]; 
							wk = 1/wk;
							if (wk<0)
								b = b+.5;
							b = b-(gw+g2)/2+gw*g2;
						}
					h2 = wh*wh;
					a2 = wk*wk;
					h4 = h2*.5;
					ex = 0;
					if (h4<150.0)
						ex = exp(-h4);
					w2 = h4*ex;
					ap = 1;
					s2 = ap-ex;
					sp = ap;
					s1 = 0;
					sn = s1;
					conex = fabs(con/wk);
					do {
						cn = ap*s2/(sn+sp);
						s1 = s1+cn;
						if (fabs(cn)<=conex)
							break;
						sn = sp;
						sp = sp+1;
						s2 = s2-w2;
						w2 = w2*h4/sp;
						ap = -ap*a2;
					} while (1);
					t = (atan(wk)-wk*s1)/twopi;
label40:
					b = b+sgn*t;
				}
				if (is>=0)
					break;
			} while(ak[0]!=0);
		}
		else if (r[0]>=0)
			if (ah[0]>=ak[0])
				b = 2*gh;
			else
				b = 2*gk;
		else if (ah[0]+ak[0]<0)
			b = 2*(gh+gk)-1;
	}
	
label50:
	if (b<0)
		b = 0;
	if (b>1)
		b = 1;

    /* free all memory allocated dinamically*/
	mxFree(buf1);
	mxFree(buf2);
	
	output[0]=b;
	return;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  double *a, *b, *rho, *y;
  int arows, acols, brows, bcols, rhorows, rhocols;
  
  /* Check for proper number of arguments. */
  if (nrhs != 3)
    mexErrMsgTxt("3 arguments in input required. See Help."); 
  
  /* The input must be a noncomplex scalar double.*/
  arows = mxGetM(prhs[0]);
  acols = mxGetN(prhs[0]);
  if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(arows == 1 && acols == 1)) {
    mexErrMsgTxt("Input argument 1 must be a noncomplex scalar double. See Help.");
  }
  brows = mxGetM(prhs[1]);
  bcols = mxGetN(prhs[1]);
  if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
      !(brows == 1 && bcols == 1)) {
    mexErrMsgTxt("Input argument 2 must be a noncomplex scalar double. See Help.");
  }
  rhorows = mxGetM(prhs[2]);
  rhocols = mxGetN(prhs[2]);
  if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      !(rhorows == 1 && rhocols == 1)) {
    mexErrMsgTxt("Input argument 3 must be a noncomplex scalar double.");
  }

  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
  
  /* Assign pointers to each input and output. */
  a = mxGetPr(prhs[0]);
  b = mxGetPr(prhs[1]);
  rho = mxGetPr(prhs[2]);
  y = mxGetPr(plhs[0]);
  
  /* Call the timestwo subroutine. */
  bivnor(y,a,b,rho);
}