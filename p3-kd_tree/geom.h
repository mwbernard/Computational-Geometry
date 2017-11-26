/* 
   Laura Toma, John Visentin 
*/


#ifndef __geom_h
#define __geom_h


typedef struct _point2d {
  int x,y; 
} point2D;


typedef struct _lineSegment2D {
    point2D p1, p2;
} lineSegment2D;


typedef struct _rect2D  {
    point2D origin;
    float width, height;
} rect2D;


//add any functions you might need to operate on these basic types

#endif
