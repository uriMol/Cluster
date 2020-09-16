/*
 * spmatUtils.h
 *
 *  Created on: 6 Aug 2020
 *      Author: Uri
 */

#include "spmat.h"
#include "modularityMax.h"
#ifndef SPMATUTILS_H_
#define SPMATUTILS_H_
#define IS_POSITIVE(X) ((X) > 0.00001)
/* returning the eigen vector of subSp */
double* getEigenVec(subSpmat *subSp, double *f, double *aVec, double *bVec, double *cVec, double *BVk);
/* dividing group g into 2 groups - g1 and g2 according to the division vector */
void divG1G2(double* division, int n, group* g, group** g1, group** g2);
/* given eigen vector returns the division vector */
double* divByEigen(double* eigenVec, int n);
/* given eigen vector - returns the eigen value */
double getEigenVal(double *eigenVec, subSpmat *sp, double *f, double *aVec, double *bVec, double *cVec, double *BVk);
/* given submatrix and divisio - computes and returns the modularity */
double getModularity(subSpmat *subSp, double *division, double *aVec, double *bVec, double *cVec, double *Bs, double *f);
/* given matrix sp - return the F vector */
double* getF(spmat *sp, group *g);
#endif /* SPMATUTILS_H_ */
