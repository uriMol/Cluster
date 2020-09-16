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
#define CHECK(expr, msg) if(!(expr)) { printf("\nERROR: %s\n ", msg); exit(2); }

typedef struct _spmat
{
	/* values - keeping the non-zero values of the matrix */
	int *values;
	/* colind - keeping the column indices of the non-zero's values */
	int *colind;
	/* rowptr - keeping for each row where its first non-zero value is from the values array */
	int *rowptr;
	/* ranks - keeping the ranks of each vertex */
	int *ranks;
	/* n - size of the matrix --> Matrix is size of (n*n) */
	int	n;
	/* 2 * edges OR sum of all ranks of original A */
	int	M;
	/*The value of '1' norm of the B matrix of this mat*/
	double shift;
} spmat;

typedef struct _subSpmat
{
	/* subValues - keeping the non-zero values of the submatrix */
	int *subValues;
	/* subColind - keeping the column indices of the non-zero's values in the sub matrix */
	int *subColind;
	/* subRanks - keeping the subRanks of the vertices - their inner submatrix ranks */
	int *subRanks;
	/* origRanks - pointer to original matrix ranks array */
	int *origRanks;
	/* g - the original vertices indices */
	int *g;
	/* n - submatrix size (n*n) */
	int	n;
	/* 2 * edges OR sum of all ranks of original A */
	int	M;
	/* subM - 2 * edges of the submatrix OR sum of all subranks */
	int	subM;
	/*The value of '1' norm of the original matrix (big enough for submatrix as well)*/
	double 	shift;
} subSpmat;

/* Creates the matrix A, allocating and adding all values */
spmat* spmat_setting(FILE *inputFile);

/* Creates a sbuMatrix from A according to indices in g */
void extractSubMatrix(spmat *sp, group *g, subSpmat *subSp);

/* Allocates a new array sparse matrix of size n with nnz non-zero elements */
spmat* spmat_allocate_array(int n, int nnz);

/* Getting the value C that we needs to shift */
double getShift(spmat *sp);

/* Frees all resources in original A matrix */
void freeSpmat(spmat *sp);

/* Frees all resources in Ag submatrix */
void freeSubSpmat(subSpmat *subsp);

#endif

