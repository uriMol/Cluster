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
	list *P, *O;
	group *g, *g1, *g2;
	spmat *sp;
	FILE *inputFile, *outputFile;
	clock_t start, end;
	double *aVec, *bVec, *cVec, *BVk;
	subSpmat *subSp;
	int i, *tmpRankPtr;


	/*SP_BUFF_SET();*/
	printf("\nIn: main, start");
	start = clock();


	inputFile = fopen(argv[1], "r");
	CHECKNEQ(inputFile, NULL, "open file error");

	/* allocating the array, setting it up with all values */
	sp = spmat_setting(inputFile);

	fclose(inputFile);

	tmpRankPtr = sp->ranks;
	for(i = 0; i < sp->n; i++){
		printf("node num %d has rank %d \n",i, *(tmpRankPtr++));
	}

	if(sp->M == 0){
		printf("\nIn:main, checking sp->M == 0 getEigenVec");
		return -1;
		/*TODO special case*/
	}
	/* allocating P and O */
	P = createP(sp->n);
	O = createO();
	createVectors(&BVk, &aVec, &bVec, &cVec, sp->n);
	subSp = createSubsp(sp);

	/* Starting the division while loop */

	while(P != NULL){
		g = P->g;
		P = P->next;
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

	freeAll(O, P, sp, subSp);

	printOutput(fopen(argv[2], "r"));
	/*TODO return s - the division*/
	CHECKEQ (argc, argc, "argc");
	end = clock();
	printf("\nIn: main, complete, took %f seconds", ((double)(end - start)/CLOCKS_PER_SEC));
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
		return;
	}

	division = divByEigen(eigenVec, subSp->n);
	division = modMaximization(subSp, division, g, aVec, bVec, cVec, BVk, f);
	Q = getModularity(subSp, division, aVec, bVec, cVec, BVk, f);
	if (!IS_POSITIVE(Q)){
		*g1 = g;
		*g2 = NULL;
		return;
	}

	divG1G2(eigenVec, subSp->n, g, g1, g2);
	if ((*g1)->len == 0 || (*g2)->len == 0){
		printf("1");
	}
	freeAfterDivision(f, division, eigenVec, g);
}

void exportData(FILE *outputFile, list *O)
{
	int numOfGroups, junk, tmpLen;
	list *OPtr;

	numOfGroups = countO(O);
	OPtr = O;
	junk = fwrite(&numOfGroups, sizeof(int), 1, outputFile);
	CHECKEQ(junk, 1, "writing numOfGroups");
	while(OPtr != NULL)
	{
		tmpLen = OPtr->g->len;
		junk = fwrite(&tmpLen, sizeof(int), 1, outputFile);
		CHECKEQ(junk, 1, "writing groupSize");
		junk = fwrite(OPtr->g->indexes, sizeof(int), tmpLen, outputFile);
		CHECKEQ(junk, tmpLen, "writing indexes");
		OPtr = OPtr->next;
	}
}


void createVectors(double **BVk, double **aVec, double **bVec, double **cVec, int n){
	*(BVk) = (double*)malloc(n*sizeof(double));
	CHECKNEQ(*BVk, NULL, "malloc BVk");
	*(aVec) = (double*)malloc(n*sizeof(double));
	CHECKNEQ(*aVec, NULL, "malloc aVec");
	*(bVec) = (double*)malloc(n*sizeof(double));
	CHECKNEQ(*bVec, NULL, "malloc bVec");
	*(cVec) = (double*)calloc(n, sizeof(double));
	CHECKNEQ(*cVec, NULL, "calloc cVec");
}

subSpmat* createSubsp(spmat *sp){
	int n;
	subSpmat *subSp;
	n = sp->n;
	subSp = (subSpmat*) malloc(sizeof(subSpmat));
	CHECKNEQ(subSp, NULL, "allocating subSp");
	subSp->subRanks = (int*) malloc(sizeof(int) * n);
	CHECKNEQ(subSp->subRanks, NULL, "allocating subSp->subRanks");
	subSp->M = sp->M;
	subSp->subValues = (int*) malloc(sizeof(int) * sp->M);
	CHECKNEQ(subSp->subValues, NULL, "allocating subSp->subValues");
	subSp->origRanks = sp->ranks;
	subSp->subColind = (int*) malloc(sizeof(int) * sp->M);
	CHECKNEQ(subSp->subColind, NULL, "allocating subSp->subColind");
	subSp->shift = sp->shift;
	return subSp;
}

