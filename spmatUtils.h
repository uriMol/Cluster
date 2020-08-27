/*
 * spmatUtils.h
 *
 *  Created on: 6 Aug 2020
 *      Author: Uri
 */

#include "spmat.h"
#ifndef SPMATUTILS_H_
#define SPMATUTILS_H_
#define IS_POSITIVE(X) ((X) > 0.00001)



double* getEigenVec(subSpmat *subSp, double *f);

void divG1G2(double* eigenVec, int n, group* g, group** g1, group** g2);

double* divByEigen(double* eigenVec, int n);

double getEigenVal(double *eigenVec, subSpmat *sp, double *f);

double getModularity(subSpmat *subSp, double *division);

double* getF(spmat *sp, group *g);

#endif /* SPMATUTILS_H_ */
