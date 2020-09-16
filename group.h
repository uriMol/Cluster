/*
 * group.h
 *
 *  Created on: 27 Aug 2020
 *      Author: Omer
 */
#ifndef GROUP_H_
#define GROUP_H_
typedef struct _group
{
	/*indexes - indices of the relevent vertices in the group*/
	int *indexes;
	/* len - length of g (Ng) */
	int	len;
} group;

/* removes the max index (in place i) and updates the length of the group */
void removeMaxIndex(group *g, int i);

#endif /* GROUP_H_ */
