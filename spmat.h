/*
 * spmat.h
 *
 *  Created on: 5 Aug 2020
 *      Author: Omer
 */
#include "group.h"
#include <stdio.h>
#ifndef _SPMAT_H
#define _SPMAT_H
#define CHECKEQ(k, expected , msg) if (k != expected) printf("\nERROR - %s " , msg)
#define CHECKNEQ(k, expected , msg) if (k == expected) printf("\nERROR - %s " , msg)


typedef struct _spmat {
	/*
	 * values - keeping the non-zero values of the matrix
	 * colind - keeping the column indices of the non-zero's values
	 * rowptr - keeping for each row where its first non-zero value is from the values array
	 */
	int *values;
	int *colind;
	int *rowptr;
	int *ranks;
	/* Matrix size (n*n) */
	int		n;
	/* 2 * edges OR sum of all ranks of original A */
	double		M;
	/*The value of '1' norm of the B matrix of this mat*/
	double 	shift;

} spmat;

typedef struct _subSpmat {
	/*
	 * values - keeping the non-zero values of the matrix
	 * colind - keeping the column indices of the non-zero's values
	 * rowptr - keeping for each row where its first non-zero value is from the values array
	 */
	int *subValues;
	int *subColind;
	int *subRanks;
	int *origRanks;
	/* g is the original indices */
	int *g;
	/* Matrix size (n*n) */
	int		n;
	/* 2 * edges OR sum of all ranks of original A */
	int 	M;
	int 	subM;
	/*The value of '1' norm of the B matrix of this mat*/
	double 	shift;

} subSpmat;


/* Creating the matrix, allocating and adding all values */
spmat* spmat_setting(FILE *inputFile);

/* Create a sbuMatrix from sp according to indices in g */
subSpmat* extractSubMatrix(spmat *sp, group *g);

/* Allocates a new arrays sparse matrix of size n with nnz non-zero elements */
spmat* spmat_allocate_array(int n, int nnz);

/* Getting the value C that we needs to shift */
double getShift(spmat *sp);

void freeSpmat(spmat *sp);

void freeSubSpmat(subSpmat *subsp);





#endif

