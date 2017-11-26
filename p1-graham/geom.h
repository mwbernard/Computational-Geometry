#ifndef __geom_h
#define __geom_h


typedef struct _point2d {
  int x,y; 
} point2D;


//a list of points 
typedef struct _pointNode pointNode; 
struct _pointNode  {
  point2D p;
  pointNode* next;
} ; 

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


/* **************************************** */
/* **************************************** */
// Geometry Functions

/* returns 2 times the signed area of triangle abc. The area is positive if c
   is to the left of ab, and negative if c is to the right of ab
 */
int signed_area2D(point2D a, point2D b, point2D c); 


/* return 1 if p,q,r collinear, and 0 otherwise */
int collinear(point2D p, point2D q, point2D r);


/* return 1 if c is  strictly left of ab; 0 otherwise */
int left (point2D a, point2D b, point2D c); 

double angle(point2D a, point2D b);

point2D lowest(point2D *p, int size);


#endif
