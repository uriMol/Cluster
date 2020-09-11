/*
 * modularityMax.c
 *
 *  Created on: 27 Aug 2020
 *      Author: Omer
 */
#include <stdio.h>
#include <stdlib.h>
#include "spmatUtils.h"
#include "modularityMax.h"
#include <time.h>
#include <math.h>

void reinitializeUnmoved(group *unmoved, int len);
void freeAfterModMax(group *unmoved, double *score, int *indices, double *improve, double *division);


double* modMaximization(subSpmat *subSp, double *division, group *g){
	int i, n, *indices, maxImproveIndex;
	group *unmoved;
	double deltaQ, *newDivision, *score, *improve;

	n = g->len;
	modInitialize(&unmoved, n, division, &newDivision, &score, &indices, &improve);

	do {
		reinitializeUnmoved(unmoved, n);
		for (i = 0; i < n; i++)
		{

			computeScoreVector2(score, subSp, newDivision);

			 /* finding and moving max deltaQ vertex */
			moveMaxVertex(score, unmoved, newDivision, indices, improve, i);
		}
		maxImproveIndex = findMaxImprove(improve, n);
		shiftUntilI(newDivision, maxImproveIndex, n, indices);
		if(maxImproveIndex == (n - 1))
		{
			deltaQ = 0;
		}
		else
		{
			deltaQ = improve[maxImproveIndex];
		}
	} while (IS_POSITIVE(deltaQ));
	freeAfterModMax(unmoved, score, indices, improve, division);
	return newDivision;
}

void modInitialize(group **unmoved, int len, double *divOrig, double **divNew, double **score, int **indices, double **improve)
{
	int i;
	double *newDivPtr, *divPtr;


	*unmoved = (group*) malloc(sizeof(group));
	CHECKNEQ(*unmoved, NULL, "allocating unmoved");
	(*unmoved)->len = len;
	(*unmoved)->indexes = (int*) malloc(sizeof(int) * len);
	CHECKNEQ((*unmoved)->indexes, NULL, "allocating unmoved->indexes");
	(*divNew) = (double*) malloc(sizeof(double) * len);
	CHECKNEQ(*divNew, NULL, "allocating divNew");
	(*score) = (double*) malloc(sizeof(double) * len);
	CHECKNEQ(*score, NULL, "allocating score");
	(*indices) = (int*) malloc(sizeof(int) * len);
	CHECKNEQ(*indices, NULL, "allocating indices");
	(*improve) = (double*) malloc(sizeof(double) * len);
	CHECKNEQ(*improve, NULL, "allocating improve");

	newDivPtr = *divNew;
	divPtr = divOrig;
	for (i = 0; i < len; i++)
	{
		*newDivPtr = *divPtr;
		newDivPtr++;
		divPtr++;
	}

}

void reinitializeUnmoved(group *unmoved, int len)
{
	int i, *indPtr;

	indPtr = unmoved->indexes;
	unmoved->len = len;
	for (i = 0; i < len; i++)
	{
		*indPtr = i;
		indPtr++;
	}
}

void computeScoreVector2(double *score, subSpmat *subSp, double* newDiv)
{
	double *scorePtr, scoreVal, *di, *dj, tmpVal, M;
	int i, j, n, *subColIndPtr, Ki, Kj, *gPtr, *tmpGPtr, Aij;

	M = (double) subSp->M;
	n = subSp->n;
	scorePtr = score;

	subColIndPtr = subSp->subColind;
	di = newDiv;
	gPtr = subSp->g;
	for(i = 0; i < n; i++)
	{
		Ki = subSp->origRanks[*gPtr];
		*di *= -1;
		scoreVal = 0;
		dj = newDiv;
		tmpGPtr = subSp->g;
		for(j = 0; j < n; j++)
		{
			Kj = subSp->origRanks[*tmpGPtr];
			tmpVal = 0;
			Aij = 0;
			if(*subColIndPtr == j)
			{
				Aij = 1;
				subColIndPtr++;
			}
			tmpVal = Aij - ((Ki * Kj) / M);
			tmpVal *= (*dj);
			scoreVal += tmpVal;
			dj++;
			tmpGPtr++;
		}
		scoreVal *= (4 * (*di));
		scoreVal += ((4 * Ki * Ki) / M);
		*scorePtr = scoreVal;
		*di *= -1;
		scorePtr++;
		di++;
		gPtr++;
	}

}



void computeScoreVector(double *score, subSpmat *subSp, double* newDiv, group *unmoved,
		double QZero, double *aVec, double *bVec, double *cVec, double *BVk, double *f)
{
	int j, *unmovedPtr;
	double *scorePtr;

	scorePtr = score;
	unmovedPtr = unmoved->indexes;

	for (j = 0; j < unmoved->len; j++)
	{
		newDiv[*unmovedPtr] *= -1;
		*scorePtr = getModularity(subSp, newDiv, aVec, bVec, cVec, BVk, f) - QZero;
		newDiv[*unmovedPtr] *= -1;
		unmovedPtr++;
		scorePtr++;
	}
}

int moveMaxVertex(double *score, group *unmoved, double *newDiv, int *indices, double *improve, int i)
{
	double max, tmpScore;
	int maxIndex, j, *unmovedPtr;

	maxIndex = 0;
	unmovedPtr = unmoved->indexes;
	max = score[*unmovedPtr];

	for (j = 0; j < unmoved->len; j++)
	{
		tmpScore = score[*unmovedPtr];
		if(tmpScore >= max)
		{
			max = tmpScore;
			maxIndex = j;
		}
		unmovedPtr++;
	}
	newDiv[unmoved->indexes[maxIndex]] *= -1;
	indices[i] = unmoved->indexes[maxIndex];

	if( i == 0)
	{
		improve[i] = max;
	}
	else
	{
		improve[i] = improve[i-1] + max;
	}
	removeMaxIndex(unmoved, maxIndex);
	return unmoved->indexes[maxIndex];
}

int findMaxImprove(double *improve, int n)
{
	double *improvePtr, max;
	int i, maxIndex;


	improvePtr = improve;
	max = *improvePtr;
	maxIndex = 0;
	for (i = 0; i < n; i++)
	{
		if(*improvePtr > max)
		{
			max = *improvePtr;
			maxIndex = i;
		}
		improvePtr++;
	}


	return maxIndex;
}

void shiftUntilI(double *newDiv, int i, int n, int *indices)
{
	int k;

	for(k = n - 1; k > i; k--)
	{
		newDiv[indices[k]] *= -1;
	}
}

void freeAfterModMax(group *unmoved, double *score, int *indices, double *improve, double *division){
	free(unmoved->indexes);
	free(unmoved);
	free(score);
	free(indices);
	free(improve);
	free(division);
}









