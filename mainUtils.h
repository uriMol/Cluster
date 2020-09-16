/*
 * mainUtils.h
 *
 *  Created on: 29 Aug 2020
 *      Author: Omer
 */

#include "spmat.h"
#include "spmatUtils.h"
#include "list.h"

#ifndef MAINUTILS_H_
#define MAINUTILS_H_

/*  This header includes all the functions that assist the main work flow, such as freeing resources or intializing structures */

/* free all given arguments */
void freeAll(list *O, list* P, spmat *sp, subSpmat *subSp,
		double *aVec, double *bVec, double *cVec, double *BVk);

/* free all given arguments (called only if there was division) */
void freeAfterDivision(double *f, double *division, double *eigenVec, group *g);

/* free all given arguments (called if there was no division) */
void freeBeforeDivision(double *f, double *eigenVec);

/* dividing the group g into 2 groups */
void divideG(spmat *sp, group *g, group **g1, group **g2, double *aVec, double *bVec,
		double *cVec, double *BVk, subSpmat *subSp);

/* exporting the division data to the output file */
void exportData(FILE *outputFile, list *O);

/* creating the vectors we will use in the program later */
void createVectors(double **BVk, double **aVec, double **bVec, double **cVec, int n);

/* after division - add the g1 and g2 droups to the correct lists (P \ O ) */
void moveGroupsToLists(group *g1,group *g2,list **P,list **O);

#endif /* MAINUTILS_H_ */
