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
#include <string.h>

int main(int argc, char* argv[]){
	spmat *sp;
	FILE *inputFile, *outputFile;
	int vartices, rank, *neighbors;

	inputFile = fopen(argv[1]);
	/* assert open success */

	/* allocating the array, setting it up with all values */
	sp = spmat_setting(inputFile);

	/* continue to create B matrix */




	return 0;
}
