#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void Inputs(char *CommandString[],int CommandNum, int *NUMVARS, int *NUMRUNS, int *TCOUNT,int *NUMSAMPS, int *NUMDATA,char *trackFile, char *paramsFile,char *mcFile,double *expec_empir,double *init_lambda, double *exit_conditions){
	int a,j,i0,i1,NUMFCNS;
	int REQ, checkArgs;
	double dtmp;
	FILE *fp,*fj;
	char dataFile[500];
	char initFile[500];

	REQ = 0;
    for(a = 1; a < CommandNum; a++) {
        if(CommandString[a][0] == '-') {
            switch(CommandString[a][1]) {
                case 'n': /* number cells */
                    *NUMVARS = atoi(CommandString[a] + 3); REQ ^= (1<<0); 
					break;
                case 'd': /* Size of saved, final MC sample */
                    *NUMDATA = atoi(CommandString[a] + 3);
									REQ ^= (1<<1); break;
                case 'r': /*number of mc runs */
                    *NUMRUNS = atoi(CommandString[a] + 3); REQ ^= (1<<2); break;
                case 'i': /*number of inner iterations */
                    *TCOUNT = atoi(CommandString[a] + 3); REQ ^= (1<<3); break;
                case 's': /*the size of the MC sample that it makes at each iteration */
                    *NUMSAMPS = atoi(CommandString[a] + 3); REQ ^= (1<<4); break;
                case 't': 
                    exit_conditions[3] = CLOCKS_PER_SEC * (double) atof(CommandString[a] + 3);
					 REQ ^= (1<<5); break;
                case 'D': /*inputs*/
                    strcpy(dataFile, CommandString[a] + 3); REQ ^= (1<<6); 
					break;
                case 'L': 
                    strcpy(trackFile, CommandString[a] + 3); REQ ^= (1<<7); break;
                case 'P': 
                    strcpy(paramsFile, CommandString[a] + 3); REQ ^= (1<<8); break;
                case 'M':
                    strcpy(mcFile, CommandString[a] + 3); REQ ^= (1<<9); break;
                case 'J':
                    strcpy(initFile, CommandString[a] + 3); REQ ^= (1<<10); break;
				case 'u':
                    exit_conditions[0] = atof(CommandString[a] + 3); REQ ^= (1<<11); break;
                case 'v':
                    exit_conditions[1] = atof(CommandString[a] + 3); REQ ^= (1<<12); break;
                case 'w':
                    exit_conditions[2]= atof(CommandString[a] + 3); REQ ^= (1<<13); break;
                default:
                    fprintf(stderr, "Bad option: %s\n", CommandString[a]); exit(-1);
            }
        }
    }
    checkArgs = 14;
    if(REQ != (1<<checkArgs) - 1) {
        fprintf(stderr, "Usage: %s\n", CommandString[0]);
        for(j = 0; j < checkArgs; j++) {
            fprintf(stderr, "%d ", (REQ >> j)&1);
        }
        fprintf(stderr, "\n");
        /****** more info on Usage ******/
        exit(-1);
    }
	
	NUMFCNS = *NUMVARS + *NUMVARS * (*NUMVARS - 1) / 2;
	
	/* for(a=0; a<4; a++){
	fprintf(stderr, "%f\n", exit_conditions[a]);
	} */
	
	if((fp = fopen(dataFile, "r")) == NULL) {
    	fprintf(stderr, "Bad file name: %s\n", dataFile);
    	exit(-1);
    }
	
    if(fscanf(fp, "%d", &i0) == EOF) {
        fprintf(stderr, "Not enough cors: %d.\n", 0);
        exit(-1);
        }
    assert(i0==*NUMVARS);
    for(j = 0; j < *NUMVARS; j++) {
        if(fscanf(fp, "%d %lf", &i0, &dtmp) == EOF) {
            fprintf(stderr, "Not enough cors: %d.\n", j);
            exit(-1);
        }
        expec_empir[j] = dtmp;
    }
    for( ; j < NUMFCNS; j++) {
        if(fscanf(fp, "%d %d %lf", &i0, &i1, &dtmp) == EOF) {
            fprintf(stderr, "Not enough cors: %d.\n", j);
            exit(-1);
        }
        expec_empir[j] = dtmp;
    }
    fclose(fp);
	
	if((fj = fopen(initFile, "r")) == NULL) {
    	fprintf(stderr, "File not found, init with independent?: %s\n", initFile);
		if(strcmp(initFile,"indep")==0){
			fprintf(stderr,"Initializing with Independent Model");
			for(j = 0; j < *NUMVARS; j++) {
				init_lambda[j] = log(expec_empir[j] / (1 - expec_empir[j]));
			}
			for( ; j < NUMFCNS; j++) {
				init_lambda[j] = 0;
			}	
		}
		else {
			fprintf(stderr,"Bad input: bye bye");
			exit(-1);
		}
    }
	else {
		if(fscanf(fj, "%d", &i0) == EOF) {
			fprintf(stderr, "Not enough initial lambdas: %d.\n", 0);
			exit(-1);
			}
			assert(i0==*NUMVARS);
		for(j = 0; j < *NUMVARS; j++) {
			if(fscanf(fj, "%d %lf", &i0, &dtmp) == EOF) {
				fprintf(stderr, "Not enough initial lambdas: %d.\n", j);
				exit(-1);
			}
			init_lambda[j] = dtmp;
		}
		for( ; j < NUMFCNS; j++) {
			if(fscanf(fj, "%d %d %lf", &i0, &i1, &dtmp) == EOF) {
				fprintf(stderr, "Not enough initial lambdas: %d.\n", j);
				exit(-1);
			}
			init_lambda[j] = dtmp;
		}
		fclose(fj);
	}
}