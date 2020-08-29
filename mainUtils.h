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

void 	printOutput(FILE *outputFile);
void 	freeAll(list *O, list* P, spmat *sp);
void freeAfterDivision(subSpmat *subsp, double *f, double *division, double *eigenVec, group *g);



#endif /* MAINUTILS_H_ */