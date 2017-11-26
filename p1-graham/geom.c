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
/* **************************************** */
/* compute the convex hull of the points in p; the points on the CH are returned as a list */
pointNode* graham_scan(point2D* p, int n) {
    // pick a point with the smallest x and y
    point2D p0 = lowest(p, n);
    
    // sort the points with the lowest and the leftmost point as a pivot
    point2D *sortedPoints;
    sortedPoints = (point2D*)malloc(sizeof(point2D));
    sortedPoints = mergeSort(p, p0, n);

    Stack *S = createStack(n);

    push(S, p0);
    push(S, sortedPoints[1]);

    for (int i = 2; i < n; i++) {
        // check if p[i] is at the left of (second(S), first(S))
        if (left(second(S), first(S), sortedPoints[i]) == 1) {
            push(S, sortedPoints[i]);
        } else if (collinear(second(S), first(S), sortedPoints[i]) == 1) {
            pop(S);
            push(S, sortedPoints[i]);
        }
        
        else {
            do {
                pop(S);

            } while (left(second(S), first(S), sortedPoints[i]) != 1 && collinear(second(S), first(S), sortedPoints[i]) != 1);
            push(S, sortedPoints[i]);
        }
    }
    // makes sure the last point added doesn't cause colinear Hull
    if (S->size >=3) {
        if (collinear(first(S), second(S), third(S)) == 1) {
            point2D top = first(S);
            pop(S);
            pop(S);
            push(S, top);
        }
    }
    
    pointNode *head = NULL;
    pointNode *curr;
    
    for (int i = 0; i < S->size; i++) {

        curr = (pointNode *)malloc(sizeof(pointNode));
        curr->p = S->elements[i];
        curr->next = head;
        head = curr;   
    }

    return head;
}
/* **************************************** */
// sorts the array of points
point2D * mergeSort(point2D *p, point2D p0, int n) {
    // base case: return when the array has less than two elements
    if (n < 2) {
        return p;
    }
    // middle index
    int mid;
    mid = n/2;
    
    // divide the array into the left and the right subarrays
    point2D *L, *R;
    L = (point2D*)malloc(mid*sizeof(point2D));
    R = (point2D*)malloc((n-mid)*sizeof(point2D));
    
    for (int i = 0; i < mid; i++) L[i] = p[i];
    for (int i = mid; i < n; i++) R[i - mid] = p[i];
    
    // recursion on left and right
    L = mergeSort(L, p0, mid);
    R = mergeSort(R, p0, n - mid);
    
    // merge left and right
    p = merge(p, L, mid, R, n - mid, p0);
    
    // free the memory allocated for L and R
    free(L);
    free(R);


    return p;
}
/* **************************************** */
// merge function to be used by Mergesort. Sorts based on angle measure of the two points.
point2D* merge(point2D *p, point2D *L, int lSize, point2D *R, int rSize, point2D p0) {
    // i = p array index, l = left array index, r = right array index
    int i, l, r;
    i = 0;
    l = 0;
    r = 0;
    
    point2D p1;
    point2D p2;
    
    double p1Angle;
    double p2Angle;
    
    double p1Dist, p2Dist;

    while (l < lSize && r < rSize) {
        p1 = L[l];
        p2 = R[r];

        p1Angle = angle(p0, p1);
        p2Angle = angle(p0, p2);
        
        if (p1Angle < p2Angle) {
            p[i++] = L[l++];
        } else if (p1Angle == p2Angle) {
            p1Dist = (p1.x - p0.x)*(p1.x - p0.x) + (p1.y + p0.y)*(p1.y + p0.y);
            p2Dist = (p2.x - p0.x)*(p2.x - p0.x) + (p2.y + p0.y)*(p2.y + p0.y);
            
            if (p1Dist < p2Dist) {
                p[i++] = L[l++];
            } else {
                p[i++] = R[r++];
            }
        }
        else {
            p[i++] = R[r++];
        }
    }
    while (l < lSize) p[i++] = L[l++];
    while (r < rSize) p[i++] = R[r++];

    return p;
}
/* **************************************** */
// returns the tan measure of two points
double angle(point2D a, point2D b) {
    double ccwAngle;
    
    if (a.y == b.y && a.x == b.x) {
        ccwAngle = -1000;
    }
    else {
        ccwAngle = atan2((double)(b.y - a.y), (double)(b.x - a.x));
        
        if (ccwAngle < 0) {
            ccwAngle = 2*PI + ccwAngle;
        }
    }
    return ccwAngle;
}

/* **************************************** */
/* **************************************** */
/* Functions to be used by the stack */

void push(Stack *S, point2D element){
    if(S->size == S->capacity){
        printf("Stack is Full\n");
    }
    else{
        S->elements[S->size++] = element;
    }
    return;
}

void pop(Stack *S){
    if(S->size==0){
        printf("Stack is Empty\n");
        return;
    }
    else{
        S->size--;
    }
    return;
}

// returns top item from stack
point2D first(Stack *S){
    if(S->size==0){
        printf("Stack is Empty\n");
        exit(0);
    }
    /* Return the topmost element */
    return S->elements[S->size-1];
}

// returns second to top item from stack
point2D second(Stack *S) {
    if (S->size < 2) {
        printf("Not enough points in Stack\n");
        exit(0);
    }
    // return second topmost element
    return S->elements[S->size - 2];
}

// returns third topmost element
point2D third(Stack *S) {
    return S->elements[S->size - 3];
}

// creates the stack
Stack * createStack(int maxElements){
    Stack *S;
    S = (Stack *)malloc(sizeof(Stack));
    /* Initialise its properties */
    S->elements = (point2D *)malloc(sizeof(point2D)*maxElements);
    S->size = 0;
    S->capacity = maxElements;
    return S;
}

/* **************************************** */
// returns the point in p with the lowest y value. If there are multiple points that have that y value,
// the point with the lowest x value is used
point2D lowest(point2D *p, int size) {
    
    int min = 10000;
    int minX = 10000;
    point2D minPoint;

    for (int i = 0; i < size; i++) {
        if (p[i].y < min) {
            min = p[i].y;
            
            minPoint = p[i];
        }
    }
    
    for (int i = 0; i < size; i++) {
        if (p[i].y == min) {
            if (p[i].x < minX) {
                minX = p[i].x;
                minPoint = p[i];
            }
        }
    }
    return minPoint;
}
