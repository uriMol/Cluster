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
	FILE *inputFile;
	double *eigenVec, eigenVal, *division, Q;


	inputFile = fopen(argv[1], "r");
	/* assert open success */

	/* allocating the array, setting it up with all values */
	sp = spmat_setting(inputFile);
	if(sp->M == 0){
		return -1;
		/*TODO special case*/
	}

	eigenVec = getEigenVec(sp);

	eigenVal = getEigenVal(eigenVec, sp);

	if (!IS_POSITIVE(eigenVal)){
		return -1;
		/*TODO special case - return undividable*/
	}

	division = divByEigen(eigenVec, sp->n);

	Q = getModularity(sp, division);

	if (!IS_POSITIVE(Q)){
		return -1;
		/*TODO special case - return undividable*/
	}

	/*TODO return s - the division*/
	CHECKEQ (argc, argc, "argc");
	return 0;
}
