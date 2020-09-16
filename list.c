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
	CHECK(P != NULL, "allocating P");
	g = (group*)malloc(sizeof(group));
	CHECK(g != NULL, "allocating g in createP");
	g->indexes = (int*) malloc(sizeof(int)*n);
	CHECK(g->indexes != NULL, "allocating g->indexes in createP");
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
	list *O = NULL;
	return O;
}


list* listAdd(list *L, group *g){
	list *newL;
	if (g->len == 0)
	{
		printf(".");
	}
	newL = (list*)malloc(sizeof(list));
	CHECK(newL != NULL, "allocating newL in listAdd");
	newL->g = g;
	newL->next = L;
	return newL;

}

int countO(list *O)
{
	int cnt;
	list *tmp;

	tmp = O;
	cnt = 0;

	while (tmp != NULL && tmp->g != NULL)
	{
		cnt++;
		tmp = tmp->next;
	}
	return cnt;
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

void freeList(list *L){
	list *lPtr;
	while(L != NULL && L->g != NULL){
		lPtr = L;
		L = L->next;
		free(lPtr->g->indexes);
		free(lPtr->g);
		free(lPtr);
	}
	if(L != NULL) free(L);
}


