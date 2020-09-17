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
#include "list.h"
#include <string.h>
#include "mainUtils.h"

int main(int argc, char* argv[])
{
	list *P, *O, *tmpList;
	group *g, *g1, *g2;
	spmat *sp;
	FILE *inputFile, *outputFile;
	double *aVec, *bVec, *cVec, *BVk;
	subSpmat *subSp;

	inputFile = fopen(argv[1], "rb");
	CHECK(inputFile != NULL, "open file error");

	/* allocating the array, setting it up with all values */
	sp = spmat_setting(inputFile);
	fclose(inputFile);
	CHECK(sp->M != 0, "invalid graph, it has no edges --> dividing by 0");

	/* allocating P and O and creating some vectors and structures what will be in use later in the program*/
	P = createP(sp->n);
	O = createO();
	createVectors(&BVk, &aVec, &bVec, &cVec, sp->n);
	subSp = createSubsp(sp);

	/* Starting the division while loop */
	while(P != NULL)
	{
		g = P->g;
		tmpList = P;
		P = P->next;
		free(tmpList);
		divideG(sp, g, &g1, &g2, aVec, bVec, cVec, BVk, subSp);
		moveGroupsToLists(g1, g2, &P, &O);
	}

	outputFile = fopen(argv[2], "wb");
	exportData(outputFile, O);
	fclose(outputFile);
	freeAll(O, P, sp, subSp, aVec, bVec, cVec, BVk);
	CHECK(argc == argc, "asserting argc is correct");
	return 0;
}
