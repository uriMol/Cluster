/*
 * group.c
 *
 *  Created on: 27 Aug 2020
 *      Author: Omer
 */
#include "group.h"
#include <stdio.h>
#include <stdlib.h>
void removeMaxIndex(group *g, int i)
{
	printf("\nIn: removeMaxIndex, start");

	g->indexes[i] = g->indexes[g->len - 1];
	g->len--;

	printf("\nIn: removeMaxIndex, complete");
}
