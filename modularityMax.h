/*
 * modularityMax.h
 *
 *  Created on: 27 Aug 2020
 *      Author: Omer
 */
#include "group.h"
#include "spmat.h"
#ifndef MODULARITYMAX_H_
#define MODULARITYMAX_H_

void modInitialize(group **unmoved, int len, double *divOrig, double **divNew, double **score, int **indices, double **improve);
void computeScoreVector(double *score, subSpmat *subSp, double* newDiv, group *unmoved, double QZero, double *aVec, double *bVec, double *cVec, double *BVk);
void moveMaxVertex(double *score, group *unmoved, double *newDiv, int *indices, double *improve, int i);
void shiftUntilI(double *newDiv, int i, int n, int *indices);
int findMaxImprove(double *improve, int n);
double* modMaximization(subSpmat *subSp,double *division, group *g, double *aVec, double *bVec, double *cVec, double *BVk);
#endif /* MODULARITYMAX_H_ */
