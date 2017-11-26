
/* 
   LT based on code from John Visentin 
*/ 

#ifndef __kdtree_h
#define __kdtree_h


#include "geom.h"
#include <stddef.h> //for size_t


typedef struct _treeNode treeNode;

struct _treeNode {
  point2D p; /* If this is a leaf node, p represents the point
		stored in this leaf.  If this is not a leaf node, p
		represents the horizontal or vertical line stored in
		this node. For a vertical line, p.y is ignored. For
		a horizontal line, p.x is ignored
	     */
  char type; /* this can be 'h' (horizontal) or 'v' (vertical), or
		'l' (leaf) depending whether the node splits with
		a horizontal line or vertical line.  Technically
		this should be an enum.
	     */
  treeNode *left, *right, *parent; /* left/below and right/above children. And for parent. */

};


typedef struct _kdtree {
  treeNode *root;  //root of the tree
  size_t count; //number of nodes in the tree
  size_t height; //height of the tree
} kdtree;

/* returns the root node of the given tree */
treeNode* kdtree_getRoot(kdtree *tree);

/* returns the point stored by the node */
point2D treeNode_getPoint(treeNode *node);

/* create a new empty tree */
kdtree* kdtree_init();

/* create a new tree representing the given array of points */
kdtree* kdtree_build(point2D *points, int n);

/* free all memory allocated for the tree, including the tree
   itself */
void kdtree_free(kdtree *tree);

/* print out information about the tree including height, number of
 nodes, and each node in an in-order traversal */
void kdtree_print(kdtree *tree);

treeNode* kdtree_build_rec(point2D *px, point2D *py, int size, int depth);

int remove_coincident_points(point2D *p, int size);

void copyPoints(point2D *p, point2D *px, point2D *py, int size);

int compyx (const void *a, const void *b);

int compxy(const void *a, const void *b);

#endif
