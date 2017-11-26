#include "geom.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#define PI 3.14159265358979323846

/* **************************************** */
/* returns the signed area of triangle abc. The area is positive if c
   is to the left of ab, and negative if c is to the right of ab
 */
int signed_area2D(point2D a, point2D b, point2D c) {
    return (b.x - a.x)*(c.y - a.y) - (c.x - a.x)*(b.y - a.y);
}
/* **************************************** */
/* return 1 if p,q,r collinear, and 0 otherwise */
int collinear(point2D p, point2D q, point2D r) {
    return signed_area2D(p,q,r) == 0;
}
/* **************************************** */
/* return 1 if c is strictly left of ab; 0 otherwise */
int left (point2D a, point2D b, point2D c) {
    return signed_area2D(a,b,c) > 0;
}
