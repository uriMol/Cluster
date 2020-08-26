/*
 * list.c
 *
 *  Created on: 26 Aug 2020
 *      Author: Uri
 */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "list.h"
#include "spmat.h"

list* createP(int n){
	int i, *gPtr;
	group *g;
	list *P;
	P = (list*) malloc(sizeof(list));
	CHECKNEQ(P, NULL, "allocating P");
	g = (group*)malloc(sizeof(group));
	CHECKNEQ(g, NULL, "allocating g in createP");
	g->indexes = (int*) malloc(sizeof(int)*n);
	CHECKNEQ(g->indexes, NULL, "allocating g->indexes in createP");
	gPtr = g->indexes;
	g->len = n;
	for(i = 0; i < n; i++){
		*gPtr = i;
		gPtr++;
	}
	P->g = g;
	P->next = NULL;
	return P;
}

list* createO(){
	list *O = (list*) malloc(sizeof(list));
	CHECKNEQ(O, NULL, "allocating O");
	O->g = NULL;
	O->next = NULL;
	return O;
}


list* listAdd(list *L, group *g){
	list *newL;

	newL = (list*)malloc(sizeof(list));
	newL->g = g;
	newL->next = L;
	return newL;

}


void printO(list* O){
	int i;
	list *tmp;
	tmp = O;
	while (tmp->g != NULL){
		printf("\n");
		for(i = 0; i < tmp->g->len; i++){
			printf("%d, ", tmp->g->indexes[i]);
		}
		tmp = tmp->next;
	}
}

