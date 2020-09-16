/*
 * mainUtils.c
 *
 *  Created on: 29 Aug 2020
 *      Author: Omer
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "mainUtils.h"

void freeAll(list *O, list* P, spmat *sp, subSpmat *subSp,
		double *aVec, double *bVec, double *cVec, double *BVk)
{
	free(aVec);
	free(bVec);
	free(cVec);
	free(BVk);
	freeList(O);
	freeList(P);
	freeSpmat(sp);
	freeSubSpmat(subSp);
}

void freeAfterDivision(double *f, double *division, double *eigenVec, group *g)
{
	freeBeforeDivision(f, eigenVec);
	free(division);
	free(g->indexes);
	free(g);
}

void freeBeforeDivision(double *f, double *eigenVec)
{
	free(f);
	free(eigenVec);
}


