#ifndef BASICTYPE_H
#define BASICTYPE_H


#include <stdio.h>


/* this is the type for an element in the priority queue; should be
 defined appropriately */
typedef struct elemType elemType;
struct elemType {
  float dist;
  int angle, x, y;
  elemType *parent;
}; 


/* functions needed in priority queue */
float getPriority(elemType x);

void printElem(elemType x);


#endif
