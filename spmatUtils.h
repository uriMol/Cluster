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



double* getEigenVec(spmat *sp);

double* divByEigen(double* eigenVec, int n);

double getEigenVal(double *eigenVec, spmat *sp);

double getModularity(spmat *sp, double *division);

#endif /* SPMATUTILS_H_ */
