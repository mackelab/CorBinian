/**
 * Gives a random shuffled order between integers 1 and numCards across the line.
 * NUMLINES: like NUMVARS
 * NUMPERLINE: numCards (number of data to be shuffled, maintaining order
 *  across lines)
 *
 * Adapted from: http://www.cs.princeton.edu/introcs/21function/Shuffle.java.html
 *
 * RUN: shuffle.x -l 40 -c 189950 -D ../../DATA/NPcelldata01.dat > ../../DATA/NPcelldata01.shuffle.dat
 *      shuffle.x -l 20 -c 200000 -D ../../DATA/fake_data.v20.d200000.b200000.out > ../../DATA/fake_data.v20.d200000.b200000.shuffle.out
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mt19937ar.h"

int numCards;
int *shuffleArray;
 
void initShuffle(int N) {
    int i, j;
    int tmpCard;
    
    numCards = N;
    shuffleArray = (int *) malloc(numCards * sizeof(int));
    
    for(i = 0; i < numCards; i++) {
        shuffleArray[i] = i;
    }
    
    for(i = 0; i < numCards; i++) {
        j = i + (int) (genrand_real2() * (numCards - i));
        tmpCard = shuffleArray[i];
        shuffleArray[i] = shuffleArray[j];
        shuffleArray[j] = tmpCard;
    }
}

int getCard(currentCard) {
    if((currentCard >= numCards) || (currentCard < 0)) {
        fprintf(stderr, "Bad card number; need 0 <= %d < %d\n", currentCard, numCards);
        exit(-1);
    }
    
    return shuffleArray[currentCard];
}

void endShuffle() {
    free(shuffleArray);
}

int main(int argc, char *argv[]) {
    int NUMLINES, NUMPERLINE, NUMLINESOUT, STARTLINEOUT;
    int **data;
    
    int a, REQ;
    int i0, i1;
    
    FILE *fp;
    char dataFile[200];
    
    /* START INITIALIZING */
    /* start reading inputs */
    REQ = 0;
    for(a = 1; a < argc; a++) {
        if(argv[a][0] == '-') {
            switch(argv[a][1]) {
                case 'l':
                    NUMLINES = atoi(argv[a] + 3); REQ ^= (1<<0); break;
                case 'c':
                    NUMPERLINE = atoi(argv[a] + 3); REQ ^= (1<<1); break;
                case 's':
                    STARTLINEOUT = atoi(argv[a] + 3); REQ ^= (1<<2); break;
                case 'v':
                    NUMLINESOUT = atoi(argv[a] + 3); REQ ^= (1<<3); break;
                case 'D':
                    strcpy(dataFile, argv[a] + 3); REQ ^= (1<<4); break;
                default:
                    fprintf(stderr, "Bad option: %s\n", argv[a]); exit(-1);
            }
        }
    }
    if(REQ != (1<<5) - 1) {
        fprintf(stderr, "Usage: %s\n", argv[0]);
        /****** more info on Usage ******/
        exit(-1);
    }
    if(STARTLINEOUT + NUMLINESOUT > NUMLINES) {
        fprintf(stderr, "Exceeding line boundary. Start %d, Finish %d, Last line %d\n", STARTLINEOUT, STARTLINEOUT + NUMLINESOUT - 1, NUMLINES - 1);
        exit(-1);
    }
    /* done reading inputs */
    
    /* start allocating space */
    data = (int **) malloc(NUMLINES * sizeof(int *));
    for(i0 = 0; i0 < NUMLINES; i0++) {
        data[i0] = (int *) malloc(NUMPERLINE * sizeof(int));
    }
    /* done allocating space */
    /* DONE INITIALIZING */
    
    /* START SHUFFLING DATA */
    /* start reading data */
    if((fp = fopen(dataFile, "r")) == NULL) {
        fprintf(stderr, "Bad file name: %s\n", dataFile);
        exit(-1);
    }
    
    for(i0 = 0; i0 < NUMLINES; i0++) {
        for(i1 = 0; i1 < NUMPERLINE; i1++) {
            if(fscanf(fp, "%d", &a) == EOF) {
                fprintf(stderr, "Not enough data.\n");
                exit(-1);
            }
            data[i0][i1] = a;
        }
    }
    /* done reading data */
    
    /* start writing shuffled data */
    initShuffle(NUMPERLINE);
    for(i0 = STARTLINEOUT; i0 < STARTLINEOUT + NUMLINESOUT; i0++) {
        for(i1 = 0; i1 < NUMPERLINE - 1; i1++) {
            fprintf(stdout, "%d ", data[i0][getCard(i1)]);
        }
        fprintf(stdout, "%d\n", data[i0][getCard(i1)]);
    }
    /* done writing shuffled data */
    /* DONE SHUFFLING DATA */
    
    /* START FREEING */
    endShuffle();
    
    for(i0 = 0; i0 < NUMLINES; i0++) {
        free(data[i0]);
    }
    free(data);
    /* DONE FREEING */
    
    return 0;
}
