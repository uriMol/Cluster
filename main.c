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
#include "SPBufferset.h"
#include "mainUtils.h"

void 	divideG(spmat *sp, group *g, group **g1, group **g2, double *aVec, double *bVec,
		double *cVec, double *BVk, subSpmat *subSp);
void 	exportData(FILE *outputFile, list *O);
void 	createVectors(double **BVk, double **aVec, double **bVec, double **cVec, int n);
subSpmat* createSubsp(spmat *sp);



int main(int argc, char* argv[]){
	list *P, *O, *tmpList;
	group *g, *g1, *g2;
	spmat *sp;
	FILE *inputFile, *outputFile;
	clock_t start, end;
	double *aVec, *bVec, *cVec, *BVk;
	subSpmat *subSp;


	/*SP_BUFF_SET();*/
	printf("In: main, start");
	start = clock();
	inputFile = fopen(argv[1], "r");
	CHECK(inputFile != NULL, "open file error");
	/* allocating the array, setting it up with all values */
	sp = spmat_setting(inputFile);

	fclose(inputFile);

	/*	if(sp->M == 0){		printf("\nIn:main, checking sp->M == 0 getEigenVec");		return -1;		TODO special case	}	*/
	CHECK(sp->M != 0, "invalid graph, it has no edges --> dividing by 0");
	/* allocating P and O */
	P = createP(sp->n);
	O = createO();
	createVectors(&BVk, &aVec, &bVec, &cVec, sp->n);
	subSp = createSubsp(sp);

	/* Starting the division while loop */

	while(P != NULL){
		g = P->g;
		tmpList = P;
		P = P->next;
		free(tmpList);
		divideG(sp, g, &g1, &g2, aVec, bVec, cVec, BVk, subSp);
		if (g2 == NULL || g2->len == 0)
		{
			O = listAdd(O, g1);
			/*TODO free if len is 0*/
		}
		else if(g1 == NULL || g1->len == 0)
		{
			O = listAdd(O, g2);
			/*TODO free if len is 0*/
		}
		else
		{
			if(g1->len == 1)
			{
				O = listAdd(O, g1);
			}
			else
			{
				P = listAdd(P, g1);
			}
			if(g2->len == 1)
			{
				O = listAdd(O, g2);
			}
			else
			{
				P = listAdd(P, g2);
			}
		}
	}
	outputFile = fopen(argv[2], "w");
	exportData(outputFile, O);
	fclose(outputFile);

	freeAll(O, P, sp, subSp, aVec, bVec, cVec, BVk);
	CHECK(argc == argc, "asserting argc is correct");
	end = clock();
	printf("\nIn: main, complete, took %f seconds\n", ((double)(end - start)/CLOCKS_PER_SEC));
	return 0;
}

void divideG(spmat *sp, group *g, group **g1, group **g2, double *aVec,
		double *bVec, double *cVec, double *BVk, subSpmat *subSp){
	double *eigenVec, eigenVal, *division, Q, *f;

	extractSubMatrix(sp, g, subSp);
	f = getF(sp, g);
	eigenVec = getEigenVec(subSp, f, aVec, bVec, cVec, BVk);
	eigenVal = getEigenVal(eigenVec, subSp, f, aVec, bVec, cVec, BVk);

	if (!IS_POSITIVE(eigenVal)){
		*g1 = g;
		*g2 = NULL;
		freeBeforeDivision(f, eigenVec);
		return;
	}

	division = divByEigen(eigenVec, subSp->n);
	division = modMaximization(subSp, division, g);
	Q = getModularity(subSp, division, aVec, bVec, cVec, BVk, f);
	if (!IS_POSITIVE(Q)){
		*g1 = g;
		*g2 = NULL;
		freeBeforeDivision(f, eigenVec);
		free(division);
		return;
	}

	divG1G2(division, subSp->n, g, g1, g2);
	freeAfterDivision(f, division, eigenVec, g);
}

void exportData(FILE *outputFile, list *O)
{
	int numOfGroups, junk, tmpLen;
	list *OPtr;

	numOfGroups = countO(O);
	OPtr = O;
	junk = fwrite(&numOfGroups, sizeof(int), 1, outputFile);
	CHECK(junk == 1, "writing numOfGroups");
	while(OPtr != NULL)
	{
		tmpLen = OPtr->g->len;
		junk = fwrite(&tmpLen, sizeof(int), 1, outputFile);
		CHECK(junk == 1, "writing groupSize");
		junk = fwrite(OPtr->g->indexes, sizeof(int), tmpLen, outputFile);
		CHECK(junk == tmpLen, "writing indexes");
		OPtr = OPtr->next;
	}
}


void createVectors(double **BVk, double **aVec, double **bVec, double **cVec, int n){
	*(BVk) = (double*)malloc(n*sizeof(double));
	CHECK(*BVk != NULL, "malloc BVk");
	*(aVec) = (double*)malloc(n*sizeof(double));
	CHECK(*aVec != NULL, "malloc aVec");
	*(bVec) = (double*)malloc(n*sizeof(double));
	CHECK(*bVec != NULL, "malloc bVec");
	*(cVec) = (double*)calloc(n, sizeof(double));
	CHECK(*cVec != NULL, "calloc cVec");
}

subSpmat* createSubsp(spmat *sp){
	int n;
	subSpmat *subSp;
	n = sp->n;
	subSp = (subSpmat*) malloc(sizeof(subSpmat));
	CHECK(subSp != NULL, "allocating subSp");
	subSp->subRanks = (int*) malloc(sizeof(int) * n);
	CHECK(subSp->subRanks != NULL, "allocating subSp->subRanks");
	subSp->M = sp->M;
	subSp->subValues = (int*) malloc(sizeof(int) * sp->M);
	CHECK(subSp->subValues != NULL, "allocating subSp->subValues");
	subSp->origRanks = sp->ranks;
	subSp->subColind = (int*) malloc(sizeof(int) * sp->M);
	CHECK(subSp->subColind != NULL, "allocating subSp->subColind");
	subSp->shift = sp->shift;
	return subSp;
}

