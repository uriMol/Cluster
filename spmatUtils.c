/*
 * spmatUtils.c
 *
 *  Created on: 6 Aug 2020
 *      Author: Uri
 */

#include <stdio.h>
#include <stdlib.h>
#include "spmatUtils.h"
#include <time.h>
#include <math.h>
#define epsilon 0.0001

double* getRandVec(int n);
void getAVmult(const double *eigenVec, const spmat *sp, double* result);
void getRanksMult(double *eigenVec, spmat *sp, double *result);
void getShiftMult(double *eigenVec, spmat *sp, double *result);
void sumAll(double *aVec, double *bVec, double *cVec, double *result, int n);
int updateEigen(double *VBk, double *eigenVec, int n);

double* getEigenVec(spmat *sp){
	int n, smallDif;
	double *eigenVec, *VBk, *aVec, *bVec, *cVec;


	n = sp->n;
	/*Initializing with random values*/
	eigenVec = getRandVec(n);
	VBk = (double*)malloc(n*sizeof(double));
	CHECKNEQ(VBk, NULL, "malloc VBk");
	aVec = (double*)malloc(n*sizeof(double));
	CHECKNEQ(aVec, NULL, "malloc aVec");
	bVec = (double*)malloc(n*sizeof(double));
	CHECKNEQ(bVec, NULL, "malloc bVec");
	cVec = (double*)malloc(n*sizeof(double));
	CHECKNEQ(cVec, NULL, "malloc cVec");

	while(smallDif != n){
		getAVmult(eigenVec, sp, aVec);
		getRanksMult(eigenVec, sp, bVec);
		getShiftMult(eigenVec, sp, cVec);
		sumAll(aVec, bVec, cVec, VBk, n);
		smallDif = updateEigen(VBk, eigenVec, n);
	}
	free(aVec);
	free(bVec);
	free(cVec);
	free(VBk);
	return eigenVec;
}

void getAVmult(const double *v, const spmat *A, double *result){
	/*initializing */
	int i, j, *colIndex, *rowPointer;
	double *value, tmp_result;

	/*setting the variables*/
	value = A->values;
	colIndex = A->colind;
	rowPointer = A->rowptr;


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

void getRanksMult(double *eigenVec, spmat *sp, double *result){
	int n, i, j, *fastRnk, *slowRnk, tmpRnk;
	double *resPtr, sum, m, *eigPtr;

	m = sp->M;
	n = sp->n;
	slowRnk = sp->ranks;
	resPtr = result;
	for(i = 0; i < n; i++){
		tmpRnk = *slowRnk;
		eigPtr = eigenVec;
		sum = 0;
		fastRnk = sp->ranks;
		for(j = 0; j < n; j++){
			sum += (tmpRnk * (*fastRnk)*(*eigPtr))/m;
			eigPtr++;
			fastRnk++;
		}
		*resPtr = sum;
		resPtr++;
		slowRnk++;
	}
}


void getShiftMult(double *eigenVec, spmat *sp, double *result){
	int n, i;
		double c, *vecPtr, *eigenPtr;
		c = sp->shift;
		n = sp->n;
		vecPtr = result;
		eigenPtr = eigenVec;
		for (i = 0; i < n; i++){
			*vecPtr = ((*eigenPtr) * c);
			vecPtr++;
			eigenPtr++;
		}
}
void sumAll(double *aVec, double *bVec, double *cVec, double *result, int n){
	int i;
	double *aPtr, *bPtr, *cPtr, *res;
	aPtr = aVec;
	bPtr = bVec;
	cPtr = cVec;
	res = result;
	for(i = 0; i < n; i++){
		*res = *aPtr - *bPtr + *cPtr;
	}
}

int updateEigen(double *VBk, double *eigenVec, int n){
	int i, cnt;
	double *vecPtr, *eigPtr;
	vecPtr = VBk;
	eigPtr = eigenVec;
	cnt = 0;
	for(i = 0; i < n; i++){
		if(fabs(*eigPtr - *vecPtr) < epsilon) cnt++;
		*eigPtr = *vecPtr;
		eigPtr++;
		vecPtr++;
	}
	return cnt;
}

double* getRandVec(int n){
	int i;
	double *v;
	srand(time(NULL));
	v = (double*)malloc(n*sizeof(double));
	for (i = 0; i < n; i++){
		v[i] = rand();
	}
	return v;
}



