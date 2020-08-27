/*
 * main.c
 *
 *  Created on: 5 Aug 2020
 *      Author: Omer
 */


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "spmat.h"
#include "spmatUtils.h"
#include "list.h"
#include <string.h>
#include "SPBufferset.h"

void divideG(spmat *sp, group *g, group **g1, group **g2);

int main(int argc, char* argv[]){
	list *P, *O;
	group *g, *g1, *g2;
	spmat *sp;
	FILE *inputFile;

	/*SP_BUFF_SET();*/

	printf("\nIn: main. Starting main");


	inputFile = fopen(argv[1], "r");
	CHECKNEQ(inputFile, NULL, "open file error");

	/* allocating the array, setting it up with all values */

	printf("\nIn:main. calling spmat_setting");

	sp = spmat_setting(inputFile);

	printf("\nIn:main. finish calling spmat_setting");

	if(sp->M == 0){
		printf("\nIn:main. checking sp->M == 0 getEigenVec");
		return -1;
		/*TODO special case*/
	}


	/*
	 * allocating P and O
	 */
	P = createP(sp->n);

	O = createO();

	/*
	 * Starting the division while loop
	 */

	while(P != NULL){
		g = P->g;
		P = P->next;
		divideG(sp, g, &g1, &g2);
		if(g1 == NULL)
		{
			O = listAdd(O, g2);
		}
		else if (g2 == NULL)
		{
			O = listAdd(O, g1);
		}
		else
		{
			if(g1->len == 1)
			{
				O = listAdd(O, g1);
			}
			else
			{
				P = listAdd(P, g1);
			}
			if(g2->len == 1)
			{
				O = listAdd(O, g2);
			}
			else
			{
				P = listAdd(P, g2);
			}
		}
		printf("\n");
		printO(O);

	}


	/*TODO return s - the division*/
	CHECKEQ (argc, argc, "argc");

	printf("\nIn: main. finish running, the divisio is ready");
	return 0;
}

void divideG(spmat *sp, group *g, group **g1, group **g2){

	subSpmat *subSp;
	double *eigenVec, eigenVal, *division, Q, *f;


	printf("\nIn:divideG, starting.");

	subSp = extractSubMatrix(sp, g);

	f = getF(sp, g);

	eigenVec = getEigenVec(subSp, f);

	eigenVal = getEigenVal(eigenVec, subSp, f);

	if (!IS_POSITIVE(eigenVal)){
		*g1 = g;
		*g2 = NULL;
		return;
	}

	division = divByEigen(eigenVec, subSp->n);

	Q = getModularity(subSp, division);

	if (!IS_POSITIVE(Q)){
		*g1 = g;
		*g2 = NULL;
		return;
	}

	divG1G2(eigenVec, subSp->n, g, g1, g2);


}



