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

/*
 * A struct to contain P - the list of groups
 */
typedef struct list{

	group *g;

	struct list* next;



} list;

list* createP(int n);
list* createO();
list* listAdd(list *L, group *g);
void printO(list* O);




#endif /* LIST_ */
