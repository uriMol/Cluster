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

void freeAll(list *O, list* P, spmat *sp, subSpmat *subSp,
		double *aVec, double *bVec, double *cVec, double *BVk)
{
	free(aVec);
	free(bVec);
	free(cVec);
	free(BVk);
	freeList(O);
	freeList(P);
	freeSpmat(sp);
	freeSubSpmat(subSp);
}

void freeAfterDivision(double *f, double *division, double *eigenVec, group *g)
{
	freeBeforeDivision(f, eigenVec);
	free(division);
	free(g->indexes);
	free(g);
}

void freeBeforeDivision(double *f, double *eigenVec)
{
	free(f);
	free(eigenVec);
}

void divideG(spmat *sp, group *g, group **g1, group **g2, double *aVec,
		double *bVec, double *cVec, double *BVk, subSpmat *subSp)
{
	double *eigenVec, eigenVal, *division, Q, *f;

	extractSubMatrix(sp, g, subSp);
	f = getF(sp, g);
	eigenVec = getEigenVec(subSp, f, aVec, bVec, cVec, BVk);
	eigenVal = getEigenVal(eigenVec, subSp, f, aVec, bVec, cVec, BVk);
	if (!IS_POSITIVE(eigenVal))
	{
		*g1 = g;
		*g2 = NULL;
		freeBeforeDivision(f, eigenVec);
		return;
	}
	division = divByEigen(eigenVec, subSp->n);
	division = modMaximization(subSp, division, g);
	Q = getModularity(subSp, division, aVec, bVec, cVec, BVk, f);
	if (!IS_POSITIVE(Q))
	{
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

void createVectors(double **BVk, double **aVec, double **bVec, double **cVec, int n)
{
	*(BVk) = (double*)malloc(n*sizeof(double));
	CHECK(*BVk != NULL, "malloc BVk");
	*(aVec) = (double*)malloc(n*sizeof(double));
	CHECK(*aVec != NULL, "malloc aVec");
	*(bVec) = (double*)malloc(n*sizeof(double));
	CHECK(*bVec != NULL, "malloc bVec");
	*(cVec) = (double*)calloc(n, sizeof(double));
	CHECK(*cVec != NULL, "calloc cVec");
}

void moveGroupsToLists(group *g1,group *g2,list *P,list *O)
{
	if (g2 == NULL || g2->len == 0)
	{
		O = listAdd(O, g1);
	}
	else if(g1 == NULL || g1->len == 0)
	{
		O = listAdd(O, g2);
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


