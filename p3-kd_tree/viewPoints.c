 /*

 What it does:
 
 - generates a set of random points in 2D, builds a kd-tree on it and
 renders it in 2D with default orthogonal projection.
 
 - program is run as:  ./viewpoints <nbpoints>
 
 - when the user presses the space bar, it regenarates the point set of same size
 
 */

#include "rtimer.h"
#include "kdtree.h"
#include "geom.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif



GLfloat red[3] = {1.0, 0.0, 0.0};
GLfloat green[3] = {0.0, 1.0, 0.0};
GLfloat blue[3] = {0.0, 0.0, 1.0};
GLfloat black[3] = {0.0, 0.0, 0.0};
GLfloat white[3] = {1.0, 1.0, 1.0};
GLfloat gray[3] = {0.5, 0.5, 0.5};
GLfloat yellow[3] = {1.0, 1.0, 0.0};
GLfloat magenta[3] = {1.0, 0.0, 1.0};
GLfloat cyan[3] = {0.0, 1.0, 1.0};




/* forward declarations of functions */
void display(void);
void keypress(unsigned char key, int x, int y);
void reset();
void mergeSort(point2D *p, int n);
void merge(point2D *p, point2D *L, int lSize, point2D *R, int rSize);
int remove_coincident_points(point2D *p, int size);
void print_dpoints(int ndp);

/* test case vars */
int init_case = 0;

int color_case = 0;

int mCount = 0;

/* global variables */
const int WINDOWSIZE = 500;
const int POINT_SIZE  = 6.0f;

//the array of n points
point2D* points = NULL;
int n;

// the kd-tree created with the points
kdtree *tree = NULL;


/* ****************************** */
/* initialize  the array of points stored in global variable points[] with random points */
void initialize_points_random() {
    
  //allocate points 
  points = malloc(n * sizeof(point2D)); 
  assert(points);
  
  int i;
  for (i=0; i<n; i++) {
    points[i].x = (int)(.1*WINDOWSIZE)/2 + rand() % ((int)(.9*WINDOWSIZE));
    points[i].y =  (int)(.1*WINDOWSIZE)/2 + rand() % ((int)(.9*WINDOWSIZE));
  }
}


/* ****************************** */
/* initialize the array of points stored in global variable points[] with vertical points */
void initialize_points_vertical() {
    
    //allocate points
    points = malloc(n * sizeof(point2D));
    assert(points);
    
    int i;
    for (i=0; i<n; i++) {
        points[i].x = WINDOWSIZE / 2;
        points[i].y = rand() % (WINDOWSIZE / 2) + (WINDOWSIZE / 4);
    }
}


/* ****************************** */

/* initialize the array of points stored in global variable points[] with horizontal points */

void initialize_points_horizontal() {
    //allocate points
    points = malloc(n * sizeof(point2D));
    assert(points);
    int i;
    for (i=0; i<n; i++) {
        points[i].y = WINDOWSIZE / 2;
        points[i].x = rand() % (WINDOWSIZE / 2) + (WINDOWSIZE / 4);
    }
}


/* ****************************** */
/* initialize the array of points stored in global variable points[] with points at a right angle */
void initialize_points_right_angle() {

    //allocate points
    
    points = malloc(n * sizeof(point2D));
    
    assert(points);
    
    points[0].x = 2 * WINDOWSIZE / 3;
    
    points[0].y = 2 * WINDOWSIZE / 3;
    
    int i;
    
    for (i=1; i < n / 2; i++) {
        
        points[i].y = 2 * WINDOWSIZE / 3;
        
        points[i].x = rand() % (WINDOWSIZE / 3) + (WINDOWSIZE / 3);
        
    }
    
    for (i = n / 2; i < n; i++) {
        
        points[i].x = 2 * WINDOWSIZE / 3;
        
        points[i].y = rand() % (WINDOWSIZE / 3) + (WINDOWSIZE / 3);
        
    }
    
}

/* ****************************** */

/* initialize the array of points stored in global variable points[] with points in a cross */

void initialize_points_cross() {
    
    
    
    //allocate points
    
    points = malloc(n * sizeof(point2D));
    
    assert(points);
    
    
    int i;
    
    for (i=0; i < n / 2; i++) {
        
        points[i].y = WINDOWSIZE / 2;
        
        points[i].x = rand() % (WINDOWSIZE / 2) + (WINDOWSIZE / 4);
        
    }
    
    for (i = n / 2; i < n; i++) {
        
        points[i].x = WINDOWSIZE / 2;
        
        points[i].y = rand() % (WINDOWSIZE / 2) + (WINDOWSIZE / 4);
        
    }
    
}


/* ****************************** */
/* print the array of points stored in global variable points[]*/
void print_points() {
    assert(points);
    int i;
    printf("points: ");
    for (i=0; i<n; i++) {
        printf("[%d,%d] ", points[i].x, points[i].y);
    }
    printf("\n");
    fflush(stdout);  //flush stdout, weird sync happens when using gl thread
}

/* ****************************** */
int main(int argc, char** argv) {
    
    // read number of points from user
    if (argc!=2) {
        printf("usage: viewPoints <nbPoints>\n");
        exit(1);
    }
    
    n = atoi(argv[1]);
    printf("you entered n=%d\n", n);
    assert(n > 0);
    
    srand(time(NULL));
    //initialize the points and build kdtree
    initialize_points_random();
    color_case = rand() % 4;
    
    tree = kdtree_init();
    
    Rtimer rt1;
    rt_start(rt1);
    tree = kdtree_build(points, n);
    rt_stop(rt1);
    char buf [1024];
    rt_sprint(buf,rt1);
    printf("time to generate kd-tree:  %s\n\n", buf);
    fflush(stdout);
    
    // print the tree
    kdtree_print(tree);
    
    
    /* initialize GLUT  */
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(WINDOWSIZE, WINDOWSIZE);
    glutInitWindowPosition(100,100);
    glutCreateWindow(argv[0]);
    
    /* register callback functions */
    glutDisplayFunc(display);
    glutKeyboardFunc(keypress);
    
    
    /* init GL */
    /* set background color black*/
    glClearColor(1, 1, 1, 1);
    
    /* circular points */
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glPointSize(POINT_SIZE);

    /* give control to event handler */
    glutMainLoop();
    return 0;
}




/* ****************************** */
/* draw a single point */
void draw_point(point2D point)
{
    glColor3fv(cyan);
    glBegin(GL_POINTS);
    glVertex2f(point.x, point.y);
    glEnd();
}


/* ****************************** */
/* draw a line between two points */
void draw_line(lineSegment2D line)
{
    GLfloat width = 2.0;
    glColor3fv(black);
    glLineWidth(width);
    glBegin(GL_LINES);
    glEnable(GL_LINE_SMOOTH);
    
    glVertex2f(line.p1.x, line.p1.y);
    glVertex2f(line.p2.x, line.p2.y);
    glEnd();
}



/* ****************************** */
/* draw the array of points stored in global variable points[] 
   each point is drawn as a small square
*/
void draw_points(){

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  //set color 
  glColor3fv(yellow);   
  
  assert(points);
  int i;
  for (i=0; i<n; i++) {
    draw_point(points[i]); 
  }
}

// fill the rectangles formed by the lines using top left and the bottom right corner
void color_rect(int x1, int y1, int x2, int y2) {
    
    GLfloat *color = malloc(3*sizeof(GLfloat));

    if (mCount % 2 == 1) {

        switch(color_case)
        {
            case 0:
                color = blue;
                break;
                
            case 1:
                color = red;
                break;
            
            case 2:
                color = white;
                break;
        
            case 3:
                color = yellow;
                break;
        }
    } else {
    
        switch(color_case)
        {
            case 0:
                color = white;
                break;
                
            case 1:
                color = yellow;
                break;
            
            case 2:
                color = blue;
                break;
        
            case 3:
                color = red;
                break;
        }
    }

    color_case = (color_case + 1) % 4;

    
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor3fv(color);
    glRectf(x1, y1, x2, y2);
    glEnd();
}

/* ****************************** */
/* recursive draw function for drawing a tree rooted at the given node
    by recursively shrinking the boundaries of the rectangle
 */
void draw_node(treeNode *node, int x1, int y1, int x2, int y2)
{
    // when the recursion gets to a leaf, the boundaries for the rectangle are set
    if (node->type == 'l') {
        color_rect(x1, y1, x2, y2);
        return;
    }

    lineSegment2D line;
    
    if (node->type == 'v') {
        int newX;
        newX = node->p.x;
        // left half
        draw_node(node->left, x1, y1, newX, y2);
        // right half
        draw_node(node->right, newX, y1, x2, y2);
        
        line.p1.x = line.p2.x = node->p.x;
        line.p1.y = y1;
        line.p2.y = y2;
    }
    else {
        int newY;
        newY = node->p.y;
        // lower half
        draw_node(node->left, x1, y1, x2, newY);
        // upper half
        draw_node(node->right, x1, newY, x2, y2);
        
        line.p1.y = line.p2.y = node->p.y;
        line.p1.x = x1;
        line.p2.x = x2;
    }
    
    draw_line(line);
    return;
    
}

/* ****************************** */
/* draw the kd-tree stored in the global variable kdTree
 */
void draw_kdtree()
{
    assert(tree);
    draw_node(tree->root, 0, 0, WINDOWSIZE, WINDOWSIZE);
    
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3fv(black);
    glRectf(0, 0, WINDOWSIZE, WINDOWSIZE);
    glEnd();
}



/* ****************************** */
void display(void) {
    
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity(); //clear the matrix
    
    
    /* the default GL window is [-1,1]x[-1,1] with the origin in the
     center the points are in the range (0,0) to (WINSIZE,WINSIZE), so
     they need to be mapped to [-1,1]x [-1,1] */
    glScalef(2.0/WINDOWSIZE, 2.0/WINDOWSIZE, 1.0);
    glTranslatef(-WINDOWSIZE/2, -WINDOWSIZE/2, 0);

    //eventually we'll want to call the function that draws the kdtree
    draw_kdtree();
   
    //for now we just draw the input points 
    //draw_points();
 
    /* execute the drawing commands */
    glFlush();
}



void reset() {
    
  //re-initialize points
    if (points) free(points);
    color_case = rand() % 4;
    mCount = 0;
    
    switch(init_case) {
        case 0:
            initialize_points_random();
            break;
        case 1:
            initialize_points_vertical();
            break;
        case 2:
            initialize_points_horizontal();
            break;
        case 3:
            initialize_points_cross();
            break;
    }
    
    init_case = (init_case + 1) % 4;
    
  //free current tree
  if (tree) kdtree_free(tree);
  
  Rtimer rt1;
  rt_start(rt1);
  tree = kdtree_build(points, n);
  rt_stop(rt1);
  char buf [1024];
  rt_sprint(buf,rt1);
  printf("time to generate kd-tree:  %s\n\n", buf);
  fflush(stdout);
  
  // print the tree
  kdtree_print(tree);
}

void movePoints() {

    int i;

    for (i = 0; i < n; i++) {
        points[i].y += (rand() % 3) - 1;
        points[i].x += (rand() % 3) - 1;
    }

    //free current tree
    if (tree) kdtree_free(tree);
   
    tree = kdtree_build(points, n);
    
}



/* ****************************** */
void keypress(unsigned char key, int x, int y) {
    switch(key)
    {
        case ' ':
            reset();
            glutPostRedisplay();
            break;

        case 'm':
            movePoints();
            glutPostRedisplay();
            mCount += 1;
            break;
            
        case 'q':
            exit(0);
            break;
    }
}


/* Handler for window re-size event. Called back when the window first appears and
 whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
    
    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);
    
    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset
    gluOrtho2D(0.0, (GLdouble) width, 0.0, (GLdouble) height); 
}


