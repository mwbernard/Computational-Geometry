/* 

Laura Toma, based on code by John Visentin

 Special cases - stop infinite recursion
 Freeing memory
 
 
*/

#include "kdtree.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>


point2D* dpoints = NULL;
int gDepth = 0;
int gCount = 0;

/* returns the root node of the given tree */
treeNode* kdtree_getRoot(kdtree *tree) {

  assert(tree); 
  return tree->root; 
}

/* returns the point stored by the node */
point2D treeNode_getPoint(treeNode *node) {

  assert(node); 
  return node->p; 
}

/* create a new empty tree */
kdtree* kdtree_init() {
  
  kdtree* tree = malloc(sizeof(kdtree));
  assert(tree); 

  tree->root = NULL; 
  tree->count = tree->height = 0; 
  return tree; 
}


/* private function to recursively print a node and its subtree in
   order
static void treeNode_print(treeNode * node) {

  if (node == NULL) return; 

  //if we are here, node must be valid

  //recursively print left child 
  treeNode_print(node->left); 

  //print node 
  switch (node->type) {
  case 'h':
    printf("HORIZONTAL: (y=%d)\n", node->p.y);
    break; 
  case 'v':
    printf("VERTICAL: (x=%d)\n", node->p.x); 
    break; 
  case 'l': 
    printf("LEAF: (p=(%d,%d))\n", node->p.x, node->p.y); 
    break;
  default: 
    printf("Improper tree node type\n"); 
    exit(1); 
  };

  //recursively print right child
  treeNode_print(node->right); 

}
*/


/* print out information about the tree including height, number of
 nodes, and each node in an in-order traversal */
void kdtree_print(kdtree *tree) {

  if (tree) {
    //printf("Nodes in order:\n\n");
    //treeNode_print(tree->root);
    printf("--- kdtree Info ---\n");
    printf("Height: %lu, Node Count: %lu\n", tree->height, tree->count);
    printf("---------------------\n");
    
    }
}

//private function to recursively free the subtree rooted at node
static void treeNode_free(treeNode* node) {
    if (node == NULL) {
        return;
    } else {
        if (node->left != NULL) {
            treeNode_free(node->left);
        }
        
        if (node->left != NULL) {
            treeNode_free(node->right);
        }
        
        free(node);
        return;
    }
}


/* free all memory allocated for the tree, including the tree
   itself */
void kdtree_free(kdtree *tree) {

  if (!tree) return; 
  treeNode_free(tree->root); 
  free(tree); 
}




/* create a new tree representing the given array of points */
kdtree* kdtree_build(point2D *points, int n) {
    kdtree *kdtree = NULL;
    kdtree = kdtree_init();
    
    // sort the points first by x and then by y
    qsort(points, n, sizeof(point2D), compxy);
    
    // get the array of distinct points
    int ndp = 0;
    ndp = remove_coincident_points(points, n);
    
    point2D *pointsByX, *pointsByY;
    pointsByX = malloc(ndp*sizeof(point2D));
    pointsByY = malloc(ndp*sizeof(point2D));
    copyPoints(dpoints, pointsByX, pointsByY, ndp);
    
    // pointsByX already has all the points sorted first by x then by y
    // --> just sort pointsByY first by y and then by x
    qsort(pointsByY, ndp, sizeof(point2D), compyx);

    int depth = 0;
    gDepth = 0;
    gCount = 0;
    treeNode *root;
    root = kdtree_build_rec(pointsByX, pointsByY, ndp, depth);
    
    kdtree->root = root;
    kdtree->count = gCount;
    kdtree->height = gDepth;

    free(pointsByX);
    free(pointsByY);
    return kdtree;
}

treeNode* kdtree_build_rec(point2D *px, point2D *py, int size, int depth) {
    // node to return
    treeNode *node = NULL;
    node = malloc(sizeof(treeNode));
    
    // base case 1: no point
    if (size == 0) {
        node = NULL;
        return node;
    }
    
    // base case 2: when there is only one point, make it a leaf node
    else if (size == 1) {
        node->p = px[0];
        node->type = 'l';
        node->left = NULL;
        node->right = NULL;
        gCount++;
        
        return node;
    }
    
    // base case 3: when there are two points
    else if (size == 2) {
        // leaves
        treeNode *lnode, *rnode;
        lnode = malloc(sizeof(treeNode));
        rnode = malloc(sizeof(treeNode));
        
        // vertical line
        if (depth % 2 == 0) {
            // the point with a smaller x is the left leaf
            lnode->p = px[0];
            lnode->type = 'l';
            lnode->left = NULL;
            lnode->right = NULL;
            
            // the point with a larger x is the right leaf
            rnode->p = px[1];
            rnode->type = 'l';
            rnode->left = NULL;
            rnode->right = NULL;
            
            node->type = 'v';
            node->p = px[0];
            node->left = lnode;
            node->right = rnode;

        }
        
        else {
            // the point with a smaller y is the left leaf
            lnode->p = py[0];
            lnode->type = 'l';
            lnode->left = NULL;
            lnode->right = NULL;

            // the point with a larger y is the right leaf
            rnode->p = py[1];
            rnode->type = 'l';
            rnode->left = NULL;
            rnode->right = NULL;
            
            node->type = 'h';
            node->p = py[0];
            node->left = lnode;
            node->right = rnode;
        }
        
        node->left->parent = node;
        node->right->parent = node;
        
        // 3 nodes are added: a node and two leaves
        gCount += 3;
        gDepth = depth + 1;
        return node;
    }
    
    else {
        // p1 and p2, each sorted by x and by y
        point2D *p1x, *p2x, *p1y, *p2y;
        
        int lSize = 0; // size of p1
        int rSize = 0; // size of p2
        
        if (depth % 2 == 0) {
            // vertical line
            node->type = 'v';
            
            int mx = 0;
            int my = 0;
            if (size % 2 == 0) { // left median if even number of points
                mx = px[size/2 - 1].x;
                my = px[size/2 - 1].y; // y value of the median
                
                node->p = px[size/2 - 1]; // assign the median x to the node
                
                lSize = size/2;
                rSize = size - lSize;
            } else {
                mx = px[size/2].x;
                my = px[size/2].y; // y value of the median
                
                node->p = px[size/2];
                
                lSize = size/2 + 1;
                rSize = size - lSize;
            }
            
            p1x = (point2D *)malloc(sizeof(point2D)*lSize);
            p2x = (point2D *)malloc(sizeof(point2D)*rSize);
            p1y = (point2D *)malloc(sizeof(point2D)*lSize);
            p2y = (point2D *)malloc(sizeof(point2D)*rSize);
            
            // p1 and p2 sorted by x --> take elements from px array
            for (int i = 0; i < lSize; i++) {
                p1x[i] = px[i];
            }
            
            int m = 0;
            for (int i = lSize; i < size; i++) {
                p2x[m] = px[i];
                m++;
            }
            
            // p1 and p2 sorted by y --> take elements from py array, which is sorted by y and then x
            int j = 0;
            int k = 0;
            for (int i = 0; i < size; i++) {
                if (py[i].x < mx) {
                    p1y[j] = py[i];
                    j++;
                }
                // if the current point in py has the same x value as the median
                else if (py[i].x == mx) {
                    // if the current point in py has the y value that's less than or equal to
                    // the y value of the median, add to p1
                    if (py[i].y <= my) {
                        p1y[j] = py[i];
                        j++;
                    } else { // otherwise, add to p2
                        p2y[k] = py[i];
                        k++;
                    }
                }
                else {
                    p2y[k] = py[i];
                    k++;
                }
            }
            
        } else {
            // horizontal line --> same as above, x and y reversed
            node->type = 'h';
            
            int mx = 0;
            int my = 0;
            
            if (size % 2 == 0) {
                mx = py[size/2 - 1].x;
                my = py[size/2 - 1].y;
                
                node->p = py[size/2 - 1];
                
                lSize = size/2;
                rSize = size - lSize;
            } else {
                mx = py[size/2].x;
                my = py[size/2].y;
                
                node->p = py[size/2];
                
                lSize = size/2 + 1;
                rSize = size - lSize;
            }
            
            p1x = malloc(sizeof(point2D)*lSize);
            p2x = malloc(sizeof(point2D)*rSize);
            p1y = malloc(sizeof(point2D)*lSize);
            p2y = malloc(sizeof(point2D)*rSize);
            
            for (int i = 0; i < lSize; i++) {
                p1y[i] = py[i];
            }
            
            int m = 0;
            for (int i = lSize; i < size; i++) {
                p2y[m] = py[i];
                m++;
            }
            
            int j = 0;
            int k = 0;
            
            for (int i = 0; i < size; i++) {
                if (px[i].y < my) {
                    p1x[j] = px[i];
                    j++;
                }
                else if (px[i].y == my) {
                    if (px[i].x <= mx) {
                        p1x[j] = px[i];
                        j++;
                    } else {
                        p2x[k] = px[i];
                        k++;
                    }
                } else {
                    p2x[k] = px[i];
                    k++;
                }
            }
        }
        
        int newDepth = 0;
        newDepth = depth + 1;
        
        node->left = kdtree_build_rec(p1x, p1y, lSize, newDepth);
        node->right = kdtree_build_rec(p2x, p2y, rSize, newDepth);
        node->left->parent = node;
        node->right->parent = node;
        
        free(p1x);
        free(p1y);
        free(p2x);
        free(p2y);
        
        gCount++;
        return node;
    }
}

// put the first instance of each element into a new array of distinct points named "dpoints"
int remove_coincident_points(point2D *p, int size) {
    int ndp = 1;
    point2D curr, prev;
    for (int i = 1; i < size; i++) {
        curr = p[i];
        prev = p[i - 1];
        
        if (curr.x != prev.x || curr.y != prev.y) {
            ndp++;
        }
    }
    
    dpoints = malloc(ndp * sizeof(point2D));
    dpoints[0] = p[0];
    int j = 1;
    for (int i = 1; i < size; i++) {
        curr = p[i];
        prev = p[i - 1];
        
        if (curr.x != prev.x || curr.y != prev.y) {
            dpoints[j] = curr;
            j++;
        }
    }
    
    return ndp;
}

// copy points into two arrays
void copyPoints(point2D *p, point2D *px, point2D *py, int size) {
    for (int i = 0; i < size; i++) {
        px[i] = p[i];
        py[i] = p[i];
    }
}

// comparison function for qsort by y and then x
int compyx (const void *a, const void *b) {
    const point2D *first = a;
    const point2D *second = b;
    
    int comp =  first->y - second->y;
    
    if (comp < 0)
        return -1;
    
    if (comp > 0)
        return 1;
    
    comp = first->x - second->x;
    
    return comp;
}

// comparator function for qsort by x and then y
int compxy(const void *a, const void *b) {
    const point2D *first = a;
    const point2D *second = b;
    
    int comp =  first->x - second->x;
    
    if (comp < 0)
        return -1;
    
    if (comp > 0)
        return 1;
    
    comp = first->y - second->y;
    
    return comp;
}

