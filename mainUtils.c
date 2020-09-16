/*
 * mainUtils.c
 *
 *  Created on: 29 Aug 2020
 *      Author: Omer
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "mainUtils.h"

void createTestGraph(int* vec, int len, FILE *outputFile)
{
	fwrite(vec, sizeof(int), len, outputFile);
}

void printOutput(FILE *outputFile)
{
	int numOfGroups, tmpLen, *tmpVert, i,j, *tmpPtr;
	fread(&numOfGroups, sizeof(int), 1, outputFile);
	printf("\nnum of groups: %d", numOfGroups);
	for (i = 0; i < numOfGroups; i++)
	{
		fread(&tmpLen, sizeof(int), 1, outputFile);
		printf("\ngroup number %d is length of: %d\n",i, tmpLen);
		tmpVert = (int*) malloc(sizeof(int) * tmpLen);
		fread(tmpVert, sizeof(int), tmpLen, outputFile);
		tmpPtr = tmpVert;
		for (j = 0; j < tmpLen; j++)
		{
			printf("%d, ", *tmpPtr);
			tmpPtr++;
		}
		printf("\n");
		free(tmpVert);
	}
}


void 	freeAll(list *O, list* P, spmat *sp, subSpmat *subSp,
		double *aVec, double *bVec, double *cVec, double *BVk){
	free(aVec);
	free(bVec);
	free(cVec);
	free(BVk);
	freeList(O);
	freeList(P);
	freeSpmat(sp);
	freeSubSpmat(subSp);
}

void freeAfterDivision(double *f, double *division, double *eigenVec, group *g){
	freeBeforeDivision(f, eigenVec);
	free(division);
	free(g->indexes);
	free(g);
}

void freeBeforeDivision(double *f, double *eigenVec){
	free(f);
	free(eigenVec);

}


