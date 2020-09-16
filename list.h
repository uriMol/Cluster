/*
 * list
 *
 *  Created on: 26 Aug 2020
 *      Author: Uri
 */
#include <stdio.h>
#include "spmat.h"
#ifndef LIST_
#define LIST_

/* This is the header of the list structure and all list methods */

/* A struct to contain the list of groups we are dividing */
typedef struct list
{
	/* g - keeps the group of vertices we are dividing */
	group *g;
	/* next - pointer to next list elelment */
	struct list* next;
} list;

/* Creates the list P with one group that contains all the vertices at first */
list* createP(int n);

/* Creates the list O and initializing it to be an empty list */
list* createO();

/* adding to list L the group g */
list* listAdd(list *L, group *g);

/* counting how many groups we have in O == how many groups we divided the graph */
int countO(list *O);

/* frees all the resources in list L */
void freeList(list *L);
#endif /* LIST_ */
