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
#include "SPBufferset.h"

int main(int argc, char* argv[]){
	spmat *sp;
	FILE *inputFile;
	double *eigenVec, eigenVal, *division, Q;

	/*SP_BUFF_SET();*/

	printf("\nIn: main. Starting main");


	inputFile = fopen(argv[1], "r");
	/* assert open success */

	/* allocating the array, setting it up with all values */

	printf("\nIn:main. calling spmat_setting");

	sp = spmat_setting(inputFile);

	printf("\nIn:main. finish calling spmat_setting");

	if(sp->M == 0){
		return -1;
		/*TODO special case*/
	}

	printf("\nIn:main. calling getEigenVec");

	eigenVec = getEigenVec(sp);

	printf("\nIn:main. finish getEigenVec");

	printf("\nIn: main. calling getEigenVal");

	eigenVal = getEigenVal(eigenVec, sp);

	printf("\nIn:main. finish getEigenVal");

	if (!IS_POSITIVE(eigenVal)){
		return -1;
		/*TODO special case - return undividable*/
	}

	printf("\nIn:main. calling divByEigen");

	division = divByEigen(eigenVec, sp->n);

	printf("\nIn:main. finish divByEigen");

	printf("\nIn:main. calling getModularity");

	Q = getModularity(sp, division);

	printf("\nIn: main. finish getModularity");

	if (!IS_POSITIVE(Q)){
		return -1;
		/*TODO special case - return undividable*/
	}

	/*TODO return s - the division*/
	CHECKEQ (argc, argc, "argc");

	printf("\nIn: main. finish running, the divisio is ready");
	return 0;
}
