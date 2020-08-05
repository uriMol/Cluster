/*
 * spmat.h
 *
 *  Created on: 5 Aug 2020
 *      Author: Omer
 */

#include <stdio.h>
#ifndef _SPMAT_H
#define _SPMAT_H
#define CHECK(k, expected , msg) if (k == expected) printf("ERROR - %s " , msg); exit(1)

typedef struct _spmat {
	/* Matrix size (n*n) */
	int		n;
	/* 2 * edges OR sum of all ranks */
	int		M;
	/* Adds row i the matrix. Called before any other call,
	 * exactly n times in order (i = 0 to n-1) */
	void	(*add_row)(struct _spmat *A, const double *row, int i, int rank);

	/* Frees all resources used by A */
	void	(*free)(struct _spmat *A);

	/* Multiplies matrix A by vector v, into result (result is pre-allocated) */
	void	(*mult)(const struct _spmat *A, const double *v, double *result);

	/* Private field for inner implementation.
	 * Should not be read or modified externally */
	void	*private;
} spmat;

/* Creating the matrix, allocating and adding all values */
spmat* spmat_setting(FILE *inputFile);

/* Allocates a new arrays sparse matrix of size n with nnz non-zero elements */
spmat* spmat_allocate_array(int n, int nnz);

#endif

