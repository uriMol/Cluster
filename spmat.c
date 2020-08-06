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
void arr_add_row(struct _spmat* , const double*, int);
void arr_free(spmat*);
void arr_mult(const spmat*, const double*, double*);



spmat* spmat_setting(FILE *inputFile){
	int k, i, tmpRank, vertices, ranks, junk, *tmpRow;
	spmat *sp;

	k = fread(&vertices, sizeof(int), 1, inputFile);
	CHECK(k, 1, "Reading File");

	/* checking sum of ranks*/
	ranks = 0;
	for (i = 0; i < vertices; i++)
	{
		k = fread(&tmpRank, sizeof(int), 1, inputFile);
		CHECK(k, 1, "Reading File");
		ranks += tmpRank;
		k = fread(&junk, sizeof(int), tmpRank, inputFile);
		CHECK(k, tmpRank, "Reading File");
	}

	/*Allocating the matrix */
	sp = spmat_allocate_array(vertices, ranks);

	/*now we read the matrix into spmat */
	rewind(inputFile);
	k = fread(&junk, sizeof(int) , 1, inputFile);
	CHECK(k, 1, "Reading File");
	tmpRow = (int*) malloc(vertices * sizeof(int));
	CHECK(tmpRow, NULL, "Allocating");
	for (i = 0; i < vertices; i++){
		k = fread(&tmpRank, sizeof(int) , 1, inputFile);
		CHECK(k, 1, "Reading File");
		k = fread(tmpRow, sizeof(int), tmpRank, inputFile);
		CHECK(k, tmpRank, "Reading File");
		sp->add_row(sp, tmpRow, i, tmpRank);
	}
	free(tmpRow);
	sp->shift = getShift(sp);

	return sp;
}

void arr_add_row(spmat *A, const int *row, int i, int rank){
	/* declaring all variables */
	int j, *colIndex, len;
	double const *index;
	double *valueIndex;
	array_helper *arr_help;


	/*initalizing all variables */
	j = 0;
	index = row;
	arr_help = (array_helper*) A->private;
	/*in case this is the first row */
	if (i == 0){
		arr_help->rowptr[0] = 0;
	}

	valueIndex = &(arr_help->values[arr_help ->rowptr[i]]);
	colIndex = &(arr_help -> colind[arr_help ->rowptr[i]]);

	/*inserting all edges*/
	while (j < rank){
		*valueIndex = 1;
		valueIndex++;
		*colIndex = *index;
		index++;
		colIndex++;
		j++;
	}

	arr_help->ranks[i] = rank;
	arr_help->rowptr[i+1] = (arr_help->rowptr[i] + rank);
}

spmat* spmat_allocate_array(int n, int nnz){
	spmat *sp;
	array_helper *arr_help;
	sp = (spmat*) malloc(sizeof(spmat));
	CHECK(sp, NULL, "Allocating");
	arr_help = (array_helper*) malloc(sizeof(array_helper));
	CHECK(arr_help, NULL, "Allocating");
	sp->private = (array_helper*) arr_help;

	/*assigning the functions */
	sp->add_row = arr_add_row;
	sp->free = arr_free;
	sp->mult = arr_mult;
	sp->n = n;
	sp->M = nnz;

	/* now defining the private field - using array_helper to keep all non-zero values */
	arr_help->values = (double*) malloc(nnz * sizeof(double));
	CHECK(arr_help->values, NULL, "Allocating");
	arr_help->colind = (int*) malloc(nnz * sizeof(int));
	CHECK(arr_help->colind, NULL, "Allocating");
	arr_help->rowptr = (int*) malloc((n+1) * sizeof(int));
	CHECK(arr_help->rowptr, NULL, "Allocating");
	arr_help->ranks = (int*) malloc(n * sizeof(int));
	CHECK(arr_help->ranks, NULL, "Allocating");
	return sp;
}

int getShift(spmat *sp){
	int i, j, rows, ranks, tmpRank, *rnkPtr, *colPtr, connected;
	double tmp, max;

	rows = sp->n;
	ranks = sp->M;
	colPtr = sp->private->colind;
	max = 0;
	for(i = 0; i < rows; i++){
		tmp = 0;
		rnkPtr = sp->private->ranks;
		tmpRank = rnkPtr[i];
		for(j = 0; j < tmpRank; j++){
			if (*colPtr == j){
				connected = 1;
				colPtr++;
			} else{
				connected = 0;
			}
			tmp += fabs(connected - ((tmpRank)*(*rnkPtr))/(double)ranks);
			rnkPtr++;
		}
		if (tmp >= max) max = tmp;
	}
	return max;
}



void arr_free(spmat *A){
	array_helper *arr_help = (array_helper*) A->private;
	free(arr_help->colind);
	free(arr_help->rowptr);
	free(arr_help->values);
	free(A->private);
	free(A);

}

void arr_mult(const spmat *A, const double *v, double *result){
	/*initializing */
	int i, j, *colIndex, *rowPointer;
	double *value, tmp_result;
	array_helper *arr_help;

	/*setting the variables*/
	arr_help = (array_helper*) A->private;
	value = arr_help->values;
	colIndex = arr_help->colind;
	rowPointer = arr_help->rowptr;


	/*looping and multiplying*/
	for (i = 0; i < A->n; i++){
		tmp_result = (double) 0;
		for (j = *rowPointer; j < *(rowPointer + 1); j++){
			tmp_result += (*value * v[*colIndex]);
			value++;
			colIndex++;
		}
		rowPointer++;
		*result = tmp_result;
		result++;
	}
}





