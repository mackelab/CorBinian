/**
 * A header file for finding expectations
 */
 
#ifndef EXERR_H
#define EXERR_H

void data_to_expec(int numvars, int numdata_s, int numdata_f,
                    int **data, double *expec);
void expec_diffs(int numvars, double *expecA, double *expecB,
                    double *avg_absdiff_expec1, double *avg_absdiff_expec2,
                    double *avg_absdiff_cov);

#endif
