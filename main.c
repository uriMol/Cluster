/*
 * main.c
 *
 *  Created on: 5 Aug 2020
 *      Author: Omer
 */


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "spmat.h"
#include "spmatUtils.h"
#include <string.h>

int main(int argc, char* argv[]){
	spmat *sp;
	FILE *inputFile, *outputFile;
	int vartices, rank, *neighbors;
	double *eigenVec;

	inputFile = fopen(argv[1]);
	/* assert open success */

	/* allocating the array, setting it up with all values */
	sp = spmat_setting(inputFile);
	if(sp->M == 0){
		return -1;
		/*TODO special case*/
	}

	eigenVec = getEigenVec(sp);




	return 0;
}
