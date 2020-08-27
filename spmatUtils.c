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
void getAVmult(const double *eigenVec, subSpmat *sp, double* result);
void getRanksMult(double *eigenVec, subSpmat *sp, double *result);
void getShiftandFMult(double *eigenVec, subSpmat *sp,double *f, double *result);
void sumAll(double *aVec, double *bVec, double *cVec, double *result, int n);
void normalize(double *vec, int n);
double vecDot(double *aVec, double *bVec, int n);
int updateEigen(double *VBk, double *eigenVec, int n);
void printA(spmat *sp);
void printB(spmat *sp);



double* getRandVec(int n){
	int i;
	double *v;
	srand(time(NULL));

	printf("\nIn:getRandVec. starting to create random vector");

	v = (double*)malloc(n*sizeof(double));
	CHECKNEQ(v, NULL, "allocating randVec v");
	for (i = 0; i < n; i++){
		v[i] = rand();
	}

	printf("\nIn:getRandVec. finish creating the random vector");
	return v;
}

double* getEigenVec(subSpmat *subSp, double* f){
	int n, smallDif;
	double *eigenVec, *VBk, *aVec, *bVec, *cVec;

	printf("\nIn:getEigenVec. starting the method");
	n = subSp->n;
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

	smallDif = 0;

	while(smallDif != n){
		getAVmult(eigenVec, subSp, aVec);
		getRanksMult(eigenVec, subSp, bVec);
		getShiftandFMult(eigenVec, subSp, f, cVec);
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


void getAVmult(const double *v, subSpmat *A, double *result){
	/*initializing */
	int i, j, *colIndex, *subRnkPtr, tmpRnk;
	double *value, tmp_result, *resPtr;

	/*setting the variables*/
	value = A->subValues;
	colIndex = A->subColind;
	subRnkPtr = A->subRanks;
	resPtr = result;

	/*looping and multiplying*/
	for (i = 0; i < A->n; i++){
		tmp_result = (double) 0;
		tmpRnk = *subRnkPtr;
		for (j = 0; j < tmpRnk; j++){
			tmp_result += (*value * v[*colIndex]);
			value++;
			colIndex++;
		}
		subRnkPtr++;
		*resPtr = tmp_result;
		resPtr++;
	}
}

void getRanksMult(double *eigenVec, subSpmat *sp, double *result){
	int n, i, j, *rnkPtr, *gPtr;
	double *resPtr, *eigPtr, m, rnkConst;

	m = sp->M;
	n = sp->n;
	rnkPtr = sp->origRanks;
	eigPtr = eigenVec;
	rnkConst = 0;
	gPtr = sp->g;
	for(j = 0; j < n; j++){
		rnkConst += (rnkPtr[*gPtr])*(*eigPtr)/m;
		eigPtr++;
		gPtr++;
	}
	resPtr = result;
	gPtr = sp->g;
	for(i = 0; i < n; i++){
		*resPtr = (rnkPtr[*gPtr]) * rnkConst;
		resPtr++;
		gPtr++;
	}
}

void getShiftandFMult(double *eigenVec, subSpmat *sp, double *f, double *result){
	int n, i;
	double *fPtr, c, *resPtr, *eigenPtr;
	c = sp->shift;
	n = sp->n;
	fPtr = f;
	resPtr = result;
	eigenPtr = eigenVec;
	for (i = 0; i < n; i++){
		*resPtr = ((*eigenPtr) * (c - *fPtr));
		resPtr++;
		eigenPtr++;
		fPtr++;
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
	if(sum == 0) return;
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
		if(fabs((*eigPtr) - (*vecPtr)) < epsilon)
			{
			cnt++;
			}
		*eigPtr = *vecPtr;
		eigPtr++;
		vecPtr++;
	}
	return cnt;
}

double getEigenVal(double *eigenVec, subSpmat *sp, double *f){
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
	getShiftandFMult(eigenVec, sp, f, cVec);
	sumAll(aVec, bVec, cVec, BVk, n);

	/*
	 * TODO check what to do if (vecDot(eigenVec, eigenVec, n) == 0)
	 */

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



void divG1G2(double* eigenVec, int n, group* g, group** g1, group** g2){
	int i, g1len, g2len, *gPtr, *g1Ptr, *g2Ptr;
	double *eigPtr;

	printf("\nIn:divG1G2. starting the method");
	eigPtr = eigenVec;
	g1len = 0;
	g2len = 0;
	*g1 = (group*)malloc(sizeof(group));
	CHECKNEQ(g1, NULL, "allocating g1 in divByEigen");
	(*g1)->indexes = (int*) malloc(sizeof(int)*n);
	CHECKNEQ((*g1)->indexes, NULL, "allocating g1->indexes in divByEigen");
	(*g2) = (group*)malloc(sizeof(group));
	CHECKNEQ(g2, NULL, "allocating g2 in divByEigen");
	(*g2)->indexes = (int*) malloc(sizeof(int)*n);
	CHECKNEQ((*g2)->indexes, NULL, "allocating g2->indexes in divByEigen");
	g1Ptr = (*g1)->indexes;
	g2Ptr = (*g2)->indexes;
	gPtr = g->indexes;
	for (i = 0; i < n; i++){
		if(IS_POSITIVE(*eigPtr)){
			*g1Ptr = *gPtr;
			g1len++;
			g1Ptr++;
		} else{
			*g2Ptr = *gPtr;
			g2len++;
			g2Ptr++;
		}
		gPtr++;
		eigPtr++;
	}
	(*g1)->len = g1len;
	(*g2)->len = g2len;
	(*g1)->indexes = (int*) realloc((*g1)->indexes, sizeof(int)*g1len);
	(*g2)->indexes = (int*) realloc((*g2)->indexes, sizeof(int)*g2len);

	printf("\nIn:divByEigen. finished division of the graph");
}

double getModularity(subSpmat *subSp, double *division){
	int n;
	double *Bs, *aVec, *bVec, *cVec, res;

	printf("\nIn:getModularity, starting the method");

	n = subSp->n;
	Bs = (double*)malloc(n*sizeof(double));
	CHECKNEQ(Bs, NULL, "malloc Bs");
	aVec = (double*)malloc(n*sizeof(double));
	CHECKNEQ(aVec, NULL, "malloc aVec");
	bVec = (double*)malloc(n*sizeof(double));
	CHECKNEQ(bVec, NULL, "malloc bVec");
	cVec = (double*)calloc(n, sizeof(double));
	CHECKNEQ(cVec, NULL, "malloc cVec");
	getAVmult(division, subSp, aVec);
	getRanksMult(division, subSp, bVec);
	sumAll(aVec, bVec, cVec, Bs, n);

	res = vecDot(division, Bs, n);

	printf("\nIn:getModularity. finished running, modularity is: %f", res);

	free(aVec);
	free(bVec);
	free(cVec);
	free(Bs);
	return res;
}

double* getF(spmat *sp, group *g){
	int i, j, len, *indPtr, *colInd, tmpRnk, tmpValSum;
	double *f, tmpRnkMult, M, rnkConst, *fPtr;

	M = sp->M;

	len = g->len;
	f = (double*) malloc(sizeof(double)*len);
	CHECKNEQ(f, NULL, "allocating f");
	fPtr = f;
	rnkConst = 0;
	for (i = 0; i < len; i++)
	{
		rnkConst += sp->ranks[g->indexes[i]];
	}
	rnkConst /= M;
	for(i = 0; i < len; i++){
		/* calculate all Aij */
		tmpValSum = 0;
		tmpRnk = sp->ranks[g->indexes[i]];
		tmpRnkMult = tmpRnk * rnkConst;
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
				tmpValSum++;
				indPtr++;
				colInd++;
				j++;
				tmpRnk--;
			}
		}
		/* calcaulate ranksMult */
		*fPtr = tmpValSum - tmpRnkMult;
		fPtr++;
	}
	return f;
}
void printA(spmat *sp){
	int i, j, n, *colPtr;
	double tmp, *valPtr;

	printf("\n");
	n = sp->n;
	colPtr = sp->colind;
	valPtr = sp->values;
	for (i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			tmp = 0;
			if(*colPtr == j){
				tmp += (*valPtr);
				colPtr++;
				valPtr++;
			}
			printf("%f, ", tmp);
		}
		printf("\n");
	}

}

void printB(spmat *sp){
	int i, j, n, *colPtr, *ranks;
	double tmp, c, *valPtr, m;

	c = sp->shift;
	n = sp->n;
	ranks = sp->ranks;
	colPtr = sp->colind;
	valPtr = sp->values;
	m = sp->M;
	printf("\n");
	for (i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			tmp = 0;
			if(*colPtr == j){
				tmp += (*valPtr);
				colPtr++;
				valPtr++;
			}
			if(i == j) tmp += c;
			tmp -= (double)(ranks[i] * ranks[j])/m;
			printf("%f, ", tmp);
		}
		printf("\n");
	}
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
		printf("i = %d, in group = %f\n", i + 1, *divPtr);
		divPtr++;
		eigPtr++;
	}

	printf("\nIn:divByEigen. finished division of the graph");

	return division;
}







