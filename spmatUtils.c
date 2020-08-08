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
void normalize(double *vec, int n);
double vecDot(double *aVec, double *bVec, int n);
int updateEigen(double *VBk, double *eigenVec, int n);


double* getRandVec(int n){
	int i;
	double *v;
	srand(time(NULL));

	printf("\nIn:getRandVec. starting to create random vector");

	v = (double*)malloc(n*sizeof(double));
	for (i = 0; i < n; i++){
		v[i] = rand();
	}

	printf("\nIn:getRandVec. finish creating the random vector");
	return v;
}

double* getEigenVec(spmat *sp){
	int n, smallDif;
	double *eigenVec, *VBk, *aVec, *bVec, *cVec;

	printf("\nIn:getEigenVec. starting the method");

	n = sp->n;
	/*Initializing with random values*/

	printf("\nIn:getEigenVec. calling getRandVec");

	eigenVec = getRandVec(n);

	printf("\nIn:getEigenVec. revieced a random vec");

	VBk = (double*)malloc(n*sizeof(double));
	CHECKNEQ(VBk, NULL, "malloc VBk");
	aVec = (double*)malloc(n*sizeof(double));
	CHECKNEQ(aVec, NULL, "malloc aVec");
	bVec = (double*)malloc(n*sizeof(double));
	CHECKNEQ(bVec, NULL, "malloc bVec");
	cVec = (double*)malloc(n*sizeof(double));
	CHECKNEQ(cVec, NULL, "malloc cVec");

	printf("\nIn:getEigenVec. starting power iterations");

	while(smallDif != n){
		getAVmult(eigenVec, sp, aVec);
		getRanksMult(eigenVec, sp, bVec);
		getShiftMult(eigenVec, sp, cVec);
		sumAll(aVec, bVec, cVec, VBk, n);
		normalize(VBk, n);
		smallDif = updateEigen(VBk, eigenVec, n);
	}

	printf("\nIn:getEigenVec. finish with power iteration");

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
	int n, i, j, *rnkPtr;
	double *resPtr, m, *eigPtr, rnkConst;

	m = sp->M;
	n = sp->n;
	rnkPtr = sp->ranks;
	eigPtr = eigenVec;
	for(j = 0; j < n; j++){
		rnkConst += (*rnkPtr)*(*eigPtr)/m;
		eigPtr++;
		rnkPtr++;
	}
	rnkPtr = sp->ranks;
	resPtr = result;
	for(i = 0; i < n; i++){
		*resPtr = (*rnkPtr) * rnkConst;
		resPtr++;
		rnkPtr++;
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
		res++;
		aPtr++;
		bPtr++;
		cPtr++;
	}
}

void normalize(double *vec, int n){
	int i;
	double sum, *vecPtr;

	sum = 0;
	vecPtr = vec;
	for(i = 0; i < n; i++){
			sum += (*vecPtr)*(*vecPtr);
			vecPtr++;
	}
	vecPtr = vec;
	sum = sqrt(sum);
	for (i = 0; i < n; i++){
		*vecPtr /= sum;
		vecPtr++;
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

double getEigenVal(double *eigenVec, spmat *sp){
	int n;
	double *aVec, *bVec, *cVec, *BVk, res;

	printf("\nIn:getEigenVal. starting the method");

	n = sp->n;

	BVk = (double*)malloc(n*sizeof(double));
	CHECKNEQ(BVk, NULL, "malloc VBk");
	aVec = (double*)malloc(n*sizeof(double));
	CHECKNEQ(aVec, NULL, "malloc aVec");
	bVec = (double*)malloc(n*sizeof(double));
	CHECKNEQ(bVec, NULL, "malloc bVec");
	cVec = (double*)malloc(n*sizeof(double));
	CHECKNEQ(cVec, NULL, "malloc cVec");

	printf("\nIn:getEigenVal. calling all mult functions and sumAll");

	getAVmult(eigenVec, sp, aVec);
	getRanksMult(eigenVec, sp, bVec);
	getShiftMult(eigenVec, sp, cVec);
	sumAll(aVec, bVec, cVec, BVk, n);

	res = ((vecDot(eigenVec, BVk, n) / vecDot(eigenVec, eigenVec, n)) - sp->shift);

	printf("\nIn:getEigenVal. finish running, eigenVal is: %f", res);

	free(aVec);
	free(bVec);
	free(cVec);
	free(BVk);
	return res;

}

double vecDot(double *aVec, double *bVec, int n){
	int i;
	double sum, *aPtr, *bPtr;

	sum = 0;
	aPtr = aVec;
	bPtr = bVec;
	for (i = 0; i < n; i++){
		sum += (*aPtr) * (*bPtr);
		aPtr++;
		bPtr++;
	}
	return sum;
}



double* divByEigen(double* eigenVec, int n){
	int i;
	double *division, *divPtr, *eigPtr;

	printf("\nIn:divByEigen. starting the method");

	division = (double*)malloc(n*sizeof(double));
	divPtr = division;
	eigPtr = eigenVec;
	for (i = 0; i < n; i++){
		(IS_POSITIVE(*eigPtr)) ? (*divPtr = 1) : (*divPtr = -1);
		divPtr++;
		eigPtr++;
	}

	printf("\nIn:divByEigen. finished division of the graph");

	return division;
}

double getModularity(spmat *sp, double *division){
	int n;
	double *Bs, *aVec, *bVec, *cVec, res;

	printf("\nIn:getModularity, starting the method");

	n = sp->n;
	Bs = (double*)malloc(n*sizeof(double));
	CHECKNEQ(Bs, NULL, "malloc Bs");
	aVec = (double*)malloc(n*sizeof(double));
	CHECKNEQ(aVec, NULL, "malloc aVec");
	bVec = (double*)malloc(n*sizeof(double));
	CHECKNEQ(bVec, NULL, "malloc bVec");
	cVec = (double*)calloc(n, sizeof(double));
	CHECKNEQ(cVec, NULL, "malloc cVec");
	getAVmult(division, sp, aVec);
	getRanksMult(division, sp, bVec);
	sumAll(aVec, bVec, cVec, Bs, n);

	res = vecDot(division, Bs, n);

	printf("\nIn:getModularity. finished running, modularity is: %f", res);

	free(aVec);
	free(bVec);
	free(cVec);
	free(Bs);
	return res;
}




