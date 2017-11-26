
#ifndef __convexhull_h
#define __convexhull_h




typedef struct Stack {
  int capacity;
  int size;
  point2D *elements;

} Stack;


/* **************************************** */
/* **************************************** */
// Scan and Sort

pointNode* graham_scan(point2D *p, int n);

point2D* mergeSort(point2D *p, point2D p0, int n);

point2D* merge(point2D *p, point2D *L, int lSize, point2D *R, int rSize, point2D p0);


/* **************************************** */
/* **************************************** */
//Stack Functions

Stack * createStack(int maxElements);

void pop(Stack *S);

point2D first(Stack *S);

point2D second(Stack *S);

point2D third(Stack *S);

void push(Stack *S, point2D element);

double angle(point2D a, point2D b);

point2D lowest(point2D *p, int size);

#endif