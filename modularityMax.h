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

/* initializing all relevent structures to perform the modularity maximization algorithm */
void modInitialize(group **unmoved, int len, double *divOrig, double **divNew, double **score, int **indices, double **improve);

/* compute the score vector in the first iteration */
void computeScoreVector(double *score, double *newDiv, subSpmat *subSp, int maxImproveIndex);

/* computes the score vector using the previous iteration */
void computeScoreVector2(double *score, subSpmat *subSp, double* newDiv);

/* moving the max vertex from the unmoved group in order to never move the same vertex twice */
int moveMaxVertex(double *score, group *unmoved, double *newDiv, int *indices, double *improve, int i);

/* shifts all relevent vertices in the new division to the other group according to the mod maximiztion psuedo code */
void shiftUntilI(double *newDiv, int i, int n, int *indices);

/* finds the vertex that shifting him to the other team contributes the most to the modularity */
int findMaxImprove(double *improve, int n);

/* performing modularity maximization  according to given psuedo code */
double* modMaximization(subSpmat *subSp, double *division, group *g);

#endif /* MODULARITYMAX_H_ */
