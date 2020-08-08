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

	printf("\nIn: spmat_setting. Starting spmat_setting");

	k = fread(&vertices, sizeof(int), 1, inputFile);
	CHECKEQ(k, 1, "Reading File");

	junk = (int*)malloc(sizeof(int) * vertices);

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

	printf("\nIn: spmat_setting. Calling spmat_allocate_array");

	sp = spmat_allocate_array(vertices, ranks);

	printf("\nIn: spmat_swtting. Finish allocating array");

	/*now we read the matrix into spmat */
	rewind(inputFile);
	k = fread(junk, sizeof(int) , 1, inputFile);
	CHECKEQ(k, 1, "Reading File");
	tmpRow = (int*) malloc(vertices * sizeof(int));
	CHECKNEQ(tmpRow, NULL, "Allocating");

	printf("\nIn: spmat_setting. Adding all rows using add_row");

	for (i = 0; i < vertices; i++){
		k = fread(&tmpRank, sizeof(int) , 1, inputFile);
		CHECKEQ(k, 1, "Reading File");
		k = fread(tmpRow, sizeof(int), tmpRank, inputFile);
		CHECKEQ(k, tmpRank, "Reading File");
		add_row(sp, tmpRow, i, tmpRank);
	}

	printf("\nIn: spmat_setting. Finish adding all rows");

	free(tmpRow);
	free(junk);


	printf("\nIn:spmat_seeting. Calling getShift");

	sp->shift = getShift(sp);

	printf("\nIn: spmat_setting. Finish getShift, finish spmat_setting");

	return sp;
}

void add_row(spmat *A, const int *row, int i, int rank){
	/* declaring all variables */
	int j, *colIndex;
	int const *index;
	double *valueIndex;

	printf("\nIn:add_row. adding row number %d", i);

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
		*colIndex = *index -1;
		index++;
		colIndex++;
		j++;
	}

	A->ranks[i] = rank;
	A->rowptr[i+1] = (A->rowptr[i] + rank);

	printf("\nIn:add_row. finish adding row number %d",i);
}

spmat* spmat_allocate_array(int n, int nnz){
	spmat *sp;

	printf("\nIn:spmat_allocate_array. starting allocation");

	sp = (spmat*) malloc(sizeof(spmat));
	CHECKNEQ(sp, NULL, "Allocating");

	sp->n = n;
	sp->M = nnz;

	/* now defining the private field - using array_helper to keep all non-zero values */
	sp->values = (double*) malloc(nnz * sizeof(double));
	CHECKNEQ(sp->values, NULL, "Allocating");
	sp->colind = (int*) malloc(nnz * sizeof(int));
	CHECKNEQ(sp->colind, NULL, "Allocating");
	sp->rowptr = (int*) malloc((n+1) * sizeof(int));
	CHECKNEQ(sp->rowptr, NULL, "Allocating");
	sp->ranks = (int*) malloc(n * sizeof(int));
	CHECKNEQ(sp->ranks, NULL, "Allocating");

	printf("\nIn: spmat_allocate_array. finish allocating");

	return sp;
}

double getShift(spmat *sp){
	int i, j, rows, ranks, tmpRank, *rnkPtr, *colPtr, connected;
	double tmp, max, result;

	printf("\nIn:getShift. start the method");

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

	printf("\nIn:getShift. finish method. the shift is: %f", max);
	return max;
}



void arr_free(spmat *A){
	free(A->colind);
	free(A->rowptr);
	free(A->values);
	free(A);

}







