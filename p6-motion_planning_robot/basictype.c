#include "basictype.h"



/* functions needed in priority queue */
float getPriority(elemType x) {
  return x.dist;
} 

void printElem(elemType x) {
  printf("[%f, %d, %d, %d] ", x.dist, x.angle, x.x, x.y);
} 
