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
	group *g;
	spmat *sp;
	subSpmat *subSp;
	FILE *inputFile;
	int i;
	double *eigenVec, eigenVal, *division, Q, *f;

	/*SP_BUFF_SET();*/

	printf("\nIn: main. Starting main");


	inputFile = fopen(argv[1], "r");
	/* assert open success */

	/* allocating the array, setting it up with all values */

	printf("\nIn:main. calling spmat_setting");

	sp = spmat_setting(inputFile);

	printf("\nIn:main. finish calling spmat_setting");

	if(sp->M == 0){
		printf("\nIn:main. checking sp->M == 0 getEigenVec");
		return -1;
		/*TODO special case*/
	}

	printf("\nIn:main. calling getEigenVec");
	/*
	P = setGroups(sp);
	g = getGroupToDivide(sp, P);
	 */
	g = (group*) malloc(sizeof(group));
	g->len = sp->n;
	g->indexes = (int*)malloc(sizeof(int) * g->len);
	for (i = 0; i < g->len; i++)
	{
		g->indexes[i] = i;
	}

	subSp = extractSubMatrix(sp, g);

	f = getF(sp, g);

	/* TODO: add setting to g */
	eigenVec = getEigenVec(subSp, f);
	free(g);

	printf("\nIn:main. finish getEigenVec");

	printf("\nIn: main. calling getEigenVal");

	eigenVal = getEigenVal(eigenVec, subSp, f);

	printf("\nIn:main. finish getEigenVal");

	if (!IS_POSITIVE(eigenVal)){
		return -1;
		/*TODO special case - return undividable*/
	}

	printf("\nIn:main. calling divByEigen");

	division = divByEigen(eigenVec, sp->n);

	printf("\nIn:main. finish divByEigen");

	printf("\nIn:main. calling getModularity");

	Q = getModularity(subSp, division);

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
