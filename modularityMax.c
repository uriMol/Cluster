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

void modInitialize(group **unmoved, int len, double *divOrig, double **divNew, double **score, int **indices, double **improve)
{
	int i, *indPtr;
	double *newDivPtr, *divPtr;

	printf("\nIn: initialize, start");

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

	indPtr = (*unmoved)->indexes;
	newDivPtr = *divNew;
	divPtr = divOrig;
	for (i = 0; i < len; i++)
	{
		*indPtr = i;
		*newDivPtr = *divPtr;
		newDivPtr++;
		divPtr++;
		indPtr++;
	}

	printf("\nIn: initialize, start");
}

void computeScoreVector(double *score, subSpmat *subSp, double* newDiv, group *unmoved, double QZero)
{
	int j, *unmovedPtr;
	double *scorePtr;

	printf("\nIn: computeScoreVector, start");

	scorePtr = score;
	unmovedPtr = unmoved->indexes;

	for (j = 0; j < unmoved->len; j++)
	{
		newDiv[*unmovedPtr] *= -1;
		*scorePtr = (getModularity(subSp, newDiv)) - QZero;
		newDiv[*unmovedPtr] *= -1;
		unmovedPtr++;
		scorePtr++;
	}

	printf("\nIn: computeScoreVector, complete");
}

void moveMaxVertex(double *score, group *unmoved, double *newDiv, int *indices, double *improve, int i)
{
	double *scorePtr, max;
	int maxIndex, j, *indicesPtr, *unmovedPtr;

	printf("\nIn: moveMaxVertex, start");

	scorePtr = score;
	max = *scorePtr;
	maxIndex = 0;
	unmovedPtr = unmoved->indexes;
	indicesPtr = indices;

	for (j = 0; j < unmoved->len; j++)
	{
		if(*scorePtr > max)
		{
			max = *scorePtr;
			maxIndex = j;
		}
		scorePtr++;
		unmovedPtr++;
	}
	newDiv[unmoved->indexes[maxIndex]] *= -1;
	*indicesPtr = unmoved->indexes[maxIndex];
	indicesPtr++;

	if( i == 0)
	{
		improve[i] = max;
	}
	else
	{
		improve[i] = improve[i-1] + max;
	}
	removeMaxIndex(unmoved, maxIndex);

	printf("\nIn: moveMaxVertex, complete");
}

int findMaxImprove(double *improve, int n)
{
	double *improvePtr, max;
	int i, maxIndex;

	printf("\nIn: findMaxImprove, start");

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

	printf("\nIn: findMaxImprove, complete");

	return maxIndex;
}

void shiftUntilI(double *newDiv, int i, int n, int *indices)
{
	int k;

	printf("\nIn: shiftUntilI, start");

	for(k = n - 1; k > i; k--)
	{
		newDiv[indices[k]] *= -1;
	}

	printf("\nIn: shiftUntilI, complete");
}










