/*
 * spmat.c
 *
 *  Created on: 5 Aug 2020
 *      Author: Omer
 */
#include "spmat.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

/* declaring all future functions */
void add_row(struct _spmat* , const int*, int, int);
void arr_free(spmat*);

spmat* spmat_setting(FILE *inputFile){
	int k, i, tmpRank, vertices, ranks, *junk, *tmpRow;
	spmat *sp;

	k = fread(&vertices, sizeof(int), 1, inputFile);
	CHECKEQ(k, 1, "Reading File");

	junk = (int*)malloc(sizeof(int) * vertices);
	CHECKNEQ(junk, NULL, "malloc junk");

	/* checking sum of ranks*/
	ranks = 0;
	for (i = 0; i < vertices; i++)
	{
		k = fread(&tmpRank, sizeof(int), 1, inputFile);
		CHECKEQ(k, 1, "Reading File");
		ranks += tmpRank;
		k = fread(junk, sizeof(int), tmpRank, inputFile);
		CHECKEQ(k, tmpRank, "Reading File");
	}
	/*Allocating the matrix */
	sp = spmat_allocate_array(vertices, ranks);
	/*now we read the matrix into spmat */
	rewind(inputFile);
	k = fread(junk, sizeof(int) , 1, inputFile);
	CHECKEQ(k, 1, "Reading File");
	tmpRow = (int*) malloc(vertices * sizeof(int));
	CHECKNEQ(tmpRow, NULL, "Allocating");
	for (i = 0; i < vertices; i++){
		k = fread(&tmpRank, sizeof(int) , 1, inputFile);
		CHECKEQ(k, 1, "Reading File");
		k = fread(tmpRow, sizeof(int), tmpRank, inputFile);
		CHECKEQ(k, tmpRank, "Reading File");
		add_row(sp, tmpRow, i, tmpRank);
	}
	free(tmpRow);
	free(junk);
	sp->shift = getShift(sp);


	return sp;
}

subSpmat* extractSubMatrix(spmat *sp, group *g)
{
	int len, i, j, subM, *colInd, tmpRnk, subTmpRnk, *subRnkPtr, *subColPtr, *indPtr, *subValPtr;
	subSpmat *subSp;

	len = g->len;
	/* allocating the arrays */
	subSp = (subSpmat*) malloc(sizeof(subSpmat));
	CHECKNEQ(subSp, NULL, "allocating subSp");
	subSp->subRanks = (int*) malloc(sizeof(int) * len);
	CHECKNEQ(subSp->subRanks, NULL, "allocating subSp->subRanks");
	subRnkPtr = subSp->subRanks;
	subSp->M = sp->M;
	subSp->subValues = (int*) malloc(sizeof(int) * sp->M);
	CHECKNEQ(subSp->subValues, NULL, "allocating subSp->subValues");
	subValPtr = subSp->subValues;
	subSp->origRanks = sp->ranks;
	subSp->subColind = (int*) malloc(sizeof(int) * sp->M);
	CHECKNEQ(subSp->subColind, NULL, "allocating subSp->subColind");
	subColPtr = subSp->subColind;
	subSp->g = g->indexes;
	subSp->shift = sp->shift;
	subSp->n = len;

	subM = 0;
	for(i = 0; i < len; i++){
		/* calculate all Aij */
		subTmpRnk = 0;
		tmpRnk = sp->ranks[g->indexes[i]];
		indPtr = g->indexes;
		colInd = &(sp->colind[sp->rowptr[indPtr[i]]]);
		j = 0;
		while (tmpRnk > 0 && j < len)
		{
			if (*colInd < *indPtr)
			{
				colInd++;
				tmpRnk--;
			}
			else if(*indPtr < *colInd)
			{
				indPtr++;
				j++;
			}
			else /* they are equal */
			{
				*subValPtr = 1;
				subValPtr++;
				*subColPtr = j;
				subColPtr++;
				subTmpRnk++;
				indPtr++;
				colInd++;
				j++;
				tmpRnk--;
			}
		}
		*subRnkPtr = subTmpRnk;
		subRnkPtr++;
		subM += subTmpRnk;
	}
	subSp->subColind = (int*)realloc(subSp->subColind, subM*sizeof(int));
	CHECKNEQ(subSp->subColind, NULL, "allocating subSp->subColind");
	subSp->subValues = (int*)realloc(subSp->subValues, subM*sizeof(int));
	CHECKNEQ(subSp->subValues, NULL, "allocating subSp->subValues");
	subSp->subM = subM;


	return subSp;
}

void add_row(spmat *A, const int *row, int i, int rank){
	/* declaring all variables */
	int j, *colIndex, *valueIndex;
	int const *index;

	/*initalizing all variables */
	j = 0;
	index = row;
	/*in case this is the first row */
	if (i == 0){
		A->rowptr[0] = 0;
	}

	valueIndex = &(A->values[A ->rowptr[i]]);
	colIndex = &(A -> colind[A ->rowptr[i]]);
	/*inserting all edges*/
	while (j < rank){
		*valueIndex = 1;
		valueIndex++;
		*colIndex = *index;
		index++;
		colIndex++;
		j++;
	}

	A->ranks[i] = rank;
	A->rowptr[i+1] = (A->rowptr[i] + rank);

}

spmat* spmat_allocate_array(int n, int nnz){
	spmat *sp;


	sp = (spmat*) malloc(sizeof(spmat));
	CHECKNEQ(sp, NULL, "Allocating");

	sp->n = n;
	sp->M = nnz;

	/* now defining the private field - using array_helper to keep all non-zero values */
	sp->values = (int*) malloc(nnz * sizeof(int));
	CHECKNEQ(sp->values, NULL, "Allocating");
	sp->colind = (int*) malloc(nnz * sizeof(int));
	CHECKNEQ(sp->colind, NULL, "Allocating");
	sp->rowptr = (int*) malloc((n+1) * sizeof(int));
	CHECKNEQ(sp->rowptr, NULL, "Allocating");
	sp->ranks = (int*) malloc(n * sizeof(int));
	CHECKNEQ(sp->ranks, NULL, "Allocating");


	return sp;
}

double getShift(spmat *sp){
	int i, j, rows, ranks, tmpRank, *rnkPtr, *colPtr, connected;
	double tmp, max, result;

	rows = sp->n;
	ranks = sp->M;
	colPtr = sp->colind;
	max = 0;
	for(i = 0; i < rows; i++){
		tmp = 0;
		rnkPtr = sp->ranks;
		tmpRank = rnkPtr[i];
		for(j = 0; j < rows; j++){
			if (*colPtr == j){
				connected = 1;
				colPtr++;
			} else{
				connected = 0;
			}
			result = fabs(connected - ((tmpRank)*(*rnkPtr))/(double)ranks);
			tmp += result;
			rnkPtr++;
		}
		if (tmp >= max) max = tmp;
	}

	return max;
}



void freeSpmat(spmat *sp){
	free(sp->colind);
	free(sp->ranks);
	free(sp->rowptr);
	free(sp->values);
	free(sp);
}

void freeSubSpmat(subSpmat *subsp){
	free(subsp->subColind);
	free(subsp->subRanks);
	free(subsp->subValues);
	free(subsp);

}








