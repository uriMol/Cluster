/*
 * spmat.h
 *
 *  Created on: 5 Aug 2020
 *      Author: Omer
 */

#include <stdio.h>
#ifndef _SPMAT_H
#define _SPMAT_H
#define CHECKEQ(k, expected , msg) if (k != expected) printf("ERROR - %s " , msg)
#define CHECKNEQ(k, expected , msg) if (k == expected) printf("ERROR - %s " , msg)


typedef struct _spmat {
	/*
	 * values - keeping the non-zero values of the matrix
	 * colind - keeping the column indices of the non-zero's values
	 * rowptr - keeping for each row where its first non-zero value is from the values array
	 */
	double *values;
	int *colind;
	int *rowptr;
	int *ranks;
	/* Matrix size (n*n) */
	int		n;
	/* 2 * edges OR sum of all ranks */
	double		M;
	/*The value of '1' norm of the B matrix of this mat*/
	double 	shift;

} spmat;

/* Creating the matrix, allocating and adding all values */
spmat* spmat_setting(FILE *inputFile);

/* Allocates a new arrays sparse matrix of size n with nnz non-zero elements */
spmat* spmat_allocate_array(int n, int nnz);

/* Getting the value C that we needs to shift */
double getShift(spmat *sp);

typedef struct _group {
	/*indexes: indices of the relevent vertices in the group*/
	int *indexes;
	/* len: length of g (Ng) */
	int		len;

} group;

#endif

