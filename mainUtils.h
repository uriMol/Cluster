/*
 * mainUtils.h
 *
 *  Created on: 29 Aug 2020
 *      Author: Omer
 */

#include "spmat.h"
#include "list.h"

#ifndef MAINUTILS_H_
#define MAINUTILS_H_

/* free all given arguments */
void freeAll(list *O, list* P, spmat *sp, subSpmat *subSp,
		double *aVec, double *bVec, double *cVec, double *BVk);

/* free all given arguments (called only if there was division) */
void freeAfterDivision(double *f, double *division, double *eigenVec, group *g);

/* free all given arguments (called if there was no division) */
void freeBeforeDivision(double *f, double *eigenVec);

#endif /* MAINUTILS_H_ */
