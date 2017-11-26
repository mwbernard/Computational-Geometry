  /* mouse2.c 

  Jeonguk Choi and Martin Bernard

  The program calculates and displays visibility graph for the set of obstacle polygons, and computes the shortest path from start to end
  When s is pressed, mouse clicks creates a polygon to add.
  When e is pressed, adding polygon stops and the polygon is added to the list of obstacles.
  When v is pressed, computes and displays the visibility graph for the current set of obstacles.
  When p is pressed, computes and displays the shortest path from start to end.
  When i is pressed, program initializes obstacles, visibility graph, and path.
  */

  #include "geom.h"
  #include "convexhull.h"


  #include <stdlib.h>
  #include <stdio.h>
  #include <math.h>
  #include <assert.h>
  #include <cmath>
  #include <limits>

  #ifdef __APPLE__
  #include <GLUT/glut.h>
  #else
  #include <GL/glut.h>
  #endif


  #include <vector> 
  #include <algorithm>
  #include <set>
  #include <utility>
  #include <iterator>

  #define PI 3.14159265358979323846

  using namespace std; 



  GLfloat red[3] = {1.0, 0.0, 0.0};
  GLfloat green[3] = {0.0, 1.0, 0.0};
  GLfloat blue[3] = {0.0, 0.0, 1.0};
  GLfloat black[3] = {0.0, 0.0, 0.0};
  GLfloat white[3] = {1.0, 1.0, 1.0};
  GLfloat gray[3] = {0.5, 0.5, 0.5};
  GLfloat yellow[3] = {1.0, 1.0, 0.0};
  GLfloat magenta[3] = {1.0, 0.0, 1.0};
  GLfloat cyan[3] = {0.0, 1.0, 1.0};

  GLint fillmode = 0;



  /* forward declarations of functions */
  void display(void);
  void keypress(unsigned char key, int x, int y);
  void mousepress(int button, int state, int x, int y);
  void timerfunc(); 

  point2D intersect(point2D a, point2D b, point2D c, point2D d);
  bool polygonInside(vector<point2D> poly1, vector<point2D> poly2);
  bool polygonIntersect(vector<point2D> poly1, vector<point2D> poly2);
  int compY (const void *a, const void *b);
  int compX (const void *a, const void *b);

  /* our coordinate system is (0,0) x (WINDOWSIZE,WINDOWSIZE) where the
     origin is the lower left corner */


  /* global variables */
  const int WINDOWSIZE = 750; 

  //the current polygon 
  vector<point2D>  poly;

  // obstacles as a vector of polygons
  vector<vector<point2D> > obstacles;
  vector<lineSegment2D> visibilityGraph;
  vector<point2D> vertices;
  vector<vector<point2D> > adjacencyList;
  vector<point2D> path;

  point2D mouse;

  double start_x = 20, start_y = 20, end_x = 700, end_y = 500;

/*
  int moveBool = 0;
  int xDir = 1;
  int yDir = 1;
*/

  //coordinates of last mouse click
  double mouse_x=-10, mouse_y=-10; 
  //initialized to a point outside the window

  //when this is 1, then clicking the mouse results in those points being stored in poly
  int poly_init_mode = 0; 

  /* 
  Usage

  void glutMouseFunc(void (*func)(int button, int state, int x, int y));

  Description

  glutMouseFunc sets the mouse callback for the current window. When a
  user presses and releases mouse buttons in the window, each press and
  each release generates a mouse callback. The button parameter is one
  of GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, or GLUT_RIGHT_BUTTON. For
  systems with only two mouse buttons, it may not be possible to
  generate GLUT_MIDDLE_BUTTON callback. For systems with a single mouse
  button, it may be possible to generate only a GLUT_LEFT_BUTTON
  callback. The state parameter is either GLUT_UP or GLUT_DOWN
  indicating whether the callback was due to a release or press
  respectively. The x and y callback parameters indicate the window
  relative coordinates when the mouse button state changed. If a
  GLUT_DOWN callback for a specific button is triggered, the program can
  assume a GLUT_UP callback for the same button will be generated
  (assuming the window still has a mouse callback registered) when the
  mouse button is released even if the mouse has moved outside the
  window.
  */

  void mousepress(int button, int state, int x, int y) {
    if(state == GLUT_DOWN) {
      mouse_x = x;
      mouse_y = y;
      //(x,y) are in wndow coordinates, where the origin is in the upper
      //left corner; our reference system has the origin in lower left
      //corner, this means we have to reflect y
      mouse_y = WINDOWSIZE - mouse_y; 
      if (poly_init_mode == 1) {
        point2D p = {mouse_x, mouse_y, -1}; 
        poly.push_back(p);
      }
    }
    
    mouse.x = (int) mouse_x;
    mouse.y = (int) mouse_y;
    glutPostRedisplay();
  }




  /* ****************************** */
  /* initialize  polygon stored in global variable poly  */
  void initialize_obstacles() {

    //clear the vector, in case something was there 
    poly.clear();


    // new obstacle polygon
    int n = 5; //size of polygon 
    double rad = 100;
    double step = 2 * M_PI / n;
    point2D p;
    for (int i=0; i<n; i++) {
      p.x = WINDOWSIZE/2 + rad * cos (i * step); 
      p.y = WINDOWSIZE/2 + rad * sin (i * step); 

      //insert the segment in the array of segments 
      poly.push_back (p);
    } //for i
    
    // add the polygon into the obstacles vector
    obstacles.push_back(poly);
  }

  void draw_circle(double x, double y, GLfloat *color){
    
    glColor3fv(color);   
    int precision = 100;
    double r = 4; 
    double theta = 0;
    glBegin(GL_POLYGON);
    for(int i = 0; i < precision; i++){
      theta = i * 2 * M_PI/precision;
      glVertex2f(x + r*cos(theta), y + r*sin(theta));
    }
    glEnd();
  }

  /* ****************************** */
  int main(int argc, char** argv) {

    initialize_obstacles();

    /* initialize GLUT  */
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(WINDOWSIZE, WINDOWSIZE);
    glutInitWindowPosition(100,100);
    glutCreateWindow(argv[0]);

    /* register callback functions */
    glutDisplayFunc(display); 
    glutKeyboardFunc(keypress);
    glutMouseFunc(mousepress);
    glutIdleFunc(timerfunc); 
    /* init GL */
    /* set background color black*/
    glClearColor(0, 0, 0, 0);   
    /* here we can enable depth testing and double buffering and so
       on */

    
    /* give control to event handler */
    glutMainLoop();
    return 0;
  }


  /* ****************************** */
  /* draw the polygon */
  void draw_polygon(vector<point2D> polygon, GLfloat *color){

    if (polygon.size() == 0) return; 

    //set color 
    glColor3fv(color);

    glLineWidth(2);  
    
    int i;
    for (i=0; i<polygon.size()-1; i++) {
      glBegin(GL_LINES);
      glVertex2f(polygon[i].x, polygon[i].y); 
      glVertex2f(polygon[i+1].x, polygon[i+1].y);
      glEnd();
    }
    //render last segment between last point and forst point 
    int last=polygon.size()-1; 
    glBegin(GL_LINES);
    glVertex2f(polygon[last].x, polygon[last].y); 
    glVertex2f(polygon[0].x, polygon[0].y);
    glEnd();
  }

  void draw_obstacles() {
    if (obstacles.size() == 0) return;

    for(int i = 0; i < obstacles.size(); i++) {
      draw_polygon(obstacles.at(i), white);
    }
  }

  void draw_VG(GLfloat *color) {
    if (visibilityGraph.size() == 0) return;

    //set color 
    glColor3fv(color);

    glLineWidth(1);  
    
    int i;
    for (i=0; i<visibilityGraph.size(); i++) {
      glBegin(GL_LINES);
      glVertex2f(visibilityGraph[i].p1.x, visibilityGraph[i].p1.y); 
      glVertex2f(visibilityGraph[i].p2.x, visibilityGraph[i].p2.y);
      glEnd();
    }
  }

  void draw_path(GLfloat *color) {
    if (path.size() == 0) return;
    
    //set color 
    glColor3fv(color);

    glLineWidth(3);  
    
    int i;
    for (i = 0; i < path.size(); i++) {
      if (i == path.size() - 1) return;
      glBegin(GL_LINES);
      glVertex2f(path.at(i).x, path.at(i).y); 
      glVertex2f(path.at(i+1).x, path.at(i+1).y);
      glEnd();
    }
  }

/************************checking if the obstacles are simple***********************************/

  // Because drawing function connects the last point to the first point, there cannot be any hole
  // Therefore, the polygon is not simple only if it has any intersection in the middle of the line
bool isPolySimple(vector<point2D> polygon) {
  point2D a, b, c, d;

  // go through each pair of edges and check if there is any intersection

  for (vector<point2D>::iterator it1 = polygon.begin(); it1 != polygon.end(); ++it1) {
    a = *it1;
    if (it1 == polygon.end() - 1) {
      b = polygon.at(0);
    } else {
      b = *(it1 + 1);
    }

    for (vector<point2D>::iterator it2 = polygon.begin(); it2 != polygon.end(); ++it2) {
      if (it1 == it2) continue;

      c = *it2;
      if (it2 == polygon.end() - 1) {
        d = polygon.at(0);
      } else {
        d = *(it2 + 1);
      }

        // check intersection between ab and cd
      point2D intersection = intersect(a, b, c, d);
      // if the two edges intersect at a point that is not a vertex, the polygon is not simple
      if (intersection.x != -1 && intersection.vertex != 1) { 
        return 0;
      }
    }
  }
  return 1;
}

// get the intersection between ab and cd
point2D intersect(point2D a, point2D b, point2D c, point2D d) {
  double s, t;
  point2D p;
  double D = (b.x - a.x)*(d.y-c.y) - (d.x-c.x)*(b.y-a.y);
  point2D *pointsSorted = new point2D[4];
  pointsSorted[0] = a;
  pointsSorted[1] = b;
  pointsSorted[2] = c;
  pointsSorted[3] = d;


  // ab and cd are parallel
  // if two lines are parallel, check the order of the four points and decide if they overlap
  if (D <= 0.001 && D >= -0.001) {
    // if ab and cd are vertical
    if (b.x == a.x) {
      // they could overlap
      if (a.x == c.x) {

        qsort(pointsSorted, 4, sizeof(point2D), compY);

        if (pointsSorted[0].x == a.x && pointsSorted[0].y == a.y) {
          if (pointsSorted[1].x == b.x && pointsSorted[1].y == b.y) {
            // b comes after a, no overlap
            p.x = -1;
            p.y = -1;
            return p;
          }
          p.x = -2;
          p.y = -2;
          p.vertex = 0;
          return p;
        }
        if (pointsSorted[0].x == b.x && pointsSorted[0].y == b.y) {
          if (pointsSorted[1].x == a.x && pointsSorted[1].y == a.y) {
            // a comes after b, no overlap
            p.x = -1;
            p.y = -1;
            return p;
          }

          p.x = -2;
          p.y = -2;
          p.vertex = 0;
          return p;
        }

        if (pointsSorted[0].x == c.x && pointsSorted[0].y == c.y) {
          if (pointsSorted[1].x == d.x && pointsSorted[1].y == d.y) {
            // d comes after c, no overlap
            p.x = -1;
            p.y = -1;
            return p;
          }

          p.x = -2;
          p.y = -2;
          p.vertex = 0;
          return p;
        }

        if (pointsSorted[0].x == d.x && pointsSorted[0].y == d.y) {
          if (pointsSorted[1].x == c.x && pointsSorted[1].y == c.y) {
            // c comes after d, no overlap
            p.x = -1;
            p.y = -1;
            return p;
          }

          p.x = -2;
          p.y = -2;
          p.vertex = 0;
          return p;
        }
      } else {
        // no overlap
        p.x = -1;
        p.y = -1;
        p.vertex = 0;
        return p;
      }
    }

    // ab and cd are parallel and not vertical
    double yIntercept1 = a.y - ((double) ((double) (b.y - a.y) / (double) (b.x - a.x))) * a.x;
    double yIntercept2 = c.y - ((double) ((double) (d.y - a.y) / (double) (d.x - c.x))) * c.x;
    // if the two lines have the same y-intercept, they could overlap
    if (abs(yIntercept1 - yIntercept2) < 0.001) {
      qsort(pointsSorted, 4, sizeof(point2D), compX);

      if (pointsSorted[0].x == a.x && pointsSorted[0].y == a.y) {
        if (pointsSorted[1].x == b.x && pointsSorted[1].y == b.y) {
          p.x = -1;
          p.y = -1;
          return p;
        }

        p.x = -2;
        p.y = -2;
        p.vertex = 0;
        return p;
      }
      if (pointsSorted[0].x == b.x && pointsSorted[0].y == b.y) {
        if (pointsSorted[1].x == a.x && pointsSorted[1].y == a.y) {
          p.x = -1;
          p.y = -1;
          return p;
        }

        p.x = -2;
        p.y = -2;
        p.vertex = 0;
        return p;
      }

      if (pointsSorted[0].x == c.x && pointsSorted[0].y == c.y) {
        if (pointsSorted[1].x == d.x && pointsSorted[1].y == d.y) {
          p.x = -1;
          p.y = -1;
          return p;
        }

        p.x = -2;
        p.y = -2;
        p.vertex = 0;
        return p;
      }


      if (pointsSorted[0].x == d.x && pointsSorted[0].y == d.y) {
        if (pointsSorted[1].x == c.x && pointsSorted[1].y == c.y) {
          p.x = -1;
          p.y = -1;
          return p;
        }

        p.x = -2;
        p.y = -2;
        p.vertex = 0;
        return p;
      }  
    } else { // parallel and different y-intercepts, no overlap
      p.x = -1;
      p.y = -1;
      return p;
    }
    
  } else { // not parallel
    s = (c.x - a.x)*(d.y - c.y) + (a.y-c.y)*(d.x-c.x);
    t = (c.x-a.x)*(b.y-a.y) + (a.y-c.y)*(b.x-a.x);
    
    if (s == 0) { // intersect at a
      p.x = a.x;
      p.y = a.y;
      if (t == 0 || t == D) {
        // intersect at vertices
        p.vertex = 1;
        return p;
      } else { // intersect inside the line
        p.vertex = 0;
        return p;
      }
    } else if (s == D) { // intersect at b
      p.x = b.x;
      p.y = b.y;
      if (t == 0 || t == D) { // intersect at vertices
        p.vertex = 1;
        return p;
      } else { // intersect inside the line
        p.vertex = 0;
        return p;
      }
    } else { // intersect inside the line
      if (t == 0) { // at c
        p.x = c.x;
        p.y = c.y;
        p.vertex = 0;
        return p;
      } else if (t == D) { // at d
        p.x = d.x;
        p.y = d.y;
        p.vertex = 0;
        return p;
      } else {
        s = s / D;
        t = t / D;
        if (s > 1.0 || s < 0.0 || t > 1.0 || t < 0.0) {
          p.x = -1;
          p.y = -1;
          return p;
        }
        p.x = a.x + s * (b.x - a.x);
        p.y = a.y + s * (b.y - a.y);
        p.vertex = 0;
        return p;
      }
    }
  }
  p.x = -2;
  p.y = -2;
  return p;
}

// ray crossing + special cases, referred to the textbook
char isInside(point2D p, vector<point2D> polygon) { 
  int i, i1;
  double x;
  int rCross = 0;
  int lCross = 0;
  bool rStrad, lStrad;
  vector<point2D> shiftedPolygon;

  for (i = 0; i < polygon.size(); i++) {
    shiftedPolygon.push_back(polygon.at(i));
    shiftedPolygon.at(i).x = shiftedPolygon.at(i).x - p.x;
    shiftedPolygon.at(i).y = shiftedPolygon.at(i).y - p.y;
  }

  for (i = 0; i < polygon.size(); i++) {
    	// if p is a vertex
    if (shiftedPolygon.at(i).x == 0 && shiftedPolygon.at(i).y == 0) return 'v'; // at vertex
    i1 = (i + shiftedPolygon.size() - 1) % shiftedPolygon.size();

    rStrad = (shiftedPolygon.at(i).y > 0) != (shiftedPolygon.at(i1).y > 0);
    lStrad = (shiftedPolygon.at(i).y < 0) != (shiftedPolygon.at(i1).y < 0);

    if (rStrad || lStrad) {
      x = (shiftedPolygon.at(i).x * shiftedPolygon.at(i1).y - shiftedPolygon.at(i1).x * shiftedPolygon.at(i).y) / (double) (shiftedPolygon.at(i1).y - shiftedPolygon.at(i).y);

      if (rStrad && x > 0) rCross++;
      if (lStrad && x < 0) lCross++;
    }
  }

  if ((rCross & 1) != (lCross & 1)) return 'e'; // on the edge

  if ((rCross & 1) == 1) return 'i'; // inside

  return 'o'; // outside
}

// 1. check if every polygon in the vector obstacles is simple
// 2. check if any polygon lies inside another polygon
// 3. check if any polygon overlaps another polygon
bool isObstaclesSimple() {
  vector<point2D> poly1, poly2;

  for (int i = 0; i < obstacles.size(); i++) {
    // check if every obstacle is simple polygon
    if (!isPolySimple(obstacles.at(i))) return 0;

    poly1 = obstacles.at(i);

    // check if any polygon intersects or lies within another polygon
    for (int j = 0; j < obstacles.size(); j++) {
      if (i == j) continue;

      poly2 = obstacles.at(j);

      if (polygonInside(poly1, poly2) || polygonInside(poly2, poly1)) return 0;
      if (polygonIntersect(poly1, poly2)) return 0;
    }
  }
  return 1;
}

// check if any edge of a polygon intersects an edge of the other polygon
bool polygonIntersect(vector<point2D> poly1, vector<point2D> poly2) {
  point2D a, b, c, d; // edge ab from poly1, edge cd from poly2
  point2D intersection;

  for (int i = 0; i < poly1.size(); i++) {
    if (i == poly1.size() - 1) {
      a = poly1.at(i);
      b = poly1.at(0);
    } else {
      a = poly1.at(i);
      b = poly1.at(i + 1);
    }

    for (int j = 0; j < poly2.size(); j++) {
      if (j == poly2.size() - 1) {
        c = poly2.at(j);
        d = poly2.at(0);
      } else {
        c = poly2.at(j);
        d = poly2.at(j + 1);
      }

      intersection = intersect(a, b, c, d);

      // if the two edges intersect anywhere within the line, return true
      if (intersection.x != -1) return 1;
    }
  }
  return 0;
}

// check if any of the vertices of poly1 lies inside poly2
bool polygonInside(vector<point2D> poly1, vector<point2D> poly2) {
  for (int i = 0; i < poly1.size(); i++) {
    if (isInside(poly1.at(i), poly2) != 'o') return 1;
  }

  return 0;
}
/****************************obstacles simple check ends**************************/

/**************computing visibility graph************************/

// returns the intersections between p-vertex and each edge
vector<point2D> findIntersections(point2D p, point2D vertex) {
  vector<point2D> intersections;
  vector<point2D> polygon;

  if (p.x == vertex.x && p.y == vertex.y) return intersections;

  for (int j = 0; j < obstacles.size(); j++) {
    polygon = obstacles.at(j);

    for (int i = 0; i < polygon.size(); i++) {
      point2D intersection;
      point2D endpoint1, endpoint2;
      if (i == polygon.size() - 1) {
        endpoint1 = polygon.at(i);
        endpoint2 = polygon.at(0);
        intersection = intersect(p, vertex, endpoint1, endpoint2);
      } else {
        endpoint1 = polygon.at(i);
        endpoint2 = polygon.at(i + 1);
        intersection = intersect(p, vertex, endpoint1, endpoint2);
      }


      if (intersection.x == -1 || intersection.x == -2) { // do not intersect
        continue;
      }

      if (((intersection.x == endpoint1.x) && (intersection.y == endpoint1.y))
        || ((intersection.x == endpoint2.x) && (intersection.y == endpoint2.y))) { // intersect at a vertex
        intersection.vertex = 1;
        intersections.push_back(intersection);
      } else {
        intersection.vertex = 0;
        intersections.push_back(intersection);
      }
    }
  }
  // if vertex is a non-obstacle, add the points as intersections
  if (vertex.x == WINDOWSIZE || vertex.x == 0 || (vertex.x == start_x && vertex.y == start_y) || (vertex.x == end_x && vertex.y == end_y)) {
    vertex.vertex = 1;
    intersections.push_back(vertex);
  }

  // delete p from intersections
  for (int i = 0; i < intersections.size(); i++) {
    if (intersections.at(i).x == p.x && intersections.at(i).y == p.y) {
      intersections.erase(intersections.begin() + i);
      i--;
    }
  }

  return intersections;
}

// get visible vertices from the point p
// if the closest intersection between p-vertex and each edge is that vertex, it is visible
vector<point2D> getVisVertices(point2D p, vector<point2D> polygon) {
  vector<point2D> visVertices;

  // check each line between p and a vertex
  for(vector<point2D>::iterator it = polygon.begin(); it != polygon.end(); ++it) {
    point2D curr = *it;
    vector<point2D> edgeIntersections;

    // get all the intersections between p-curr and every edge of the obstacles
    edgeIntersections = findIntersections(p, curr);

    // check the closest intersection point from p
    int minDist = 999999999;
    point2D minDistIntersection;
    for (int i = 0; i < edgeIntersections.size(); i++) {
      int dist;
      dist = (p.x - edgeIntersections.at(i).x) * (p.x - edgeIntersections.at(i).x) +
      (p.y - edgeIntersections.at(i).y) * (p.y - edgeIntersections.at(i).y);

      if (dist < minDist) {
        minDist = dist;
        minDistIntersection = edgeIntersections.at(i);
      }
    }

    if (minDistIntersection.vertex == 1 && abs(minDistIntersection.x - curr.x) < 4 && abs(minDistIntersection.y - curr.y) < 4) {
      visVertices.push_back(minDistIntersection);
    }
  }

  return visVertices;
}

// graham scan function takes an array; copy a vector into an array
point2D *copyVec(vector<point2D> polygon, int size) {
  point2D *pointArray;
  pointArray = (point2D*)malloc(size*sizeof(point2D));

  for (int i = 0; i < size; i++) {
    pointArray[i] = polygon.at(i);
  }

  return pointArray;
}

// check if the convex hull edge ab intersects other polygons
int chIntersect(point2D a, point2D b, int index) {
  for (int i = 0; i < obstacles.size(); i++) {
    // skip the polygon that the convex hull edge belongs to
    if (i == index) continue;
    for (int j = 0; j < obstacles.at(i).size(); j++) {
      point2D c, d;

      if (j == obstacles.at(i).size() - 1) {
        c = obstacles.at(i).at(j);
        d = obstacles.at(i).at(0);
      } else {
        c = obstacles.at(i).at(j);
        d = obstacles.at(i).at(j + 1);
      }

      if (intersect(a, b, c, d).x != -1) return 1;
    }
  }
  return 0;
}

// compute visibility graph and add the visible edges as lineSegment2D into the global vector
void computeVG() {
  visibilityGraph.clear();
  path.clear();

  if (obstacles.size() < 1) return;

  if (!isObstaclesSimple()) {
    printf("The obstacles are not simple polygons. Press i to initialize the scene and s to draw new obstacles.\n");
    return;
  }

  point2D start, end;
  start.x = start_x;
  start.y = start_y;
  end.x = end_x;
  end.y = end_y;

  // check if the start and end points are inside any of the polygons
  for (int i = 0; i < obstacles.size(); i++) {
    if (isInside(start, obstacles.at(i)) != 'o' || isInside(end, obstacles.at(i)) != 'o') {
      printf("The start or the goal is inside an obstacle. Press i to initialize the scene and s to draw new obstacles.\n");
      return;
    }
  }

  // go through each point of the obstacles and compute visible vertices
  lineSegment2D visPair;
  vector<point2D> visVertices;

  for (int i = 0; i < obstacles.size(); i++) { // for each obstacle polygon
    /**** compute convex hull and walk through it to add visible edges ****/

    // turn vector into array for graham_scan function
    point2D *pointArray;
    pointArray = copyVec(obstacles.at(i), obstacles.at(i).size());

    /******* have to check if the convex hull edges intersects any of the polygon *******/
    pointNode *head = graham_scan(pointArray, obstacles.at(i).size());
    pointNode *curr = head;

    int chEdgeIntersect;
    while (curr->next != NULL) {
      chEdgeIntersect = chIntersect(curr->p, curr->next->p, i);

      if (chEdgeIntersect == 0) {
        visPair.p1 = curr->p;
        visPair.p2 = curr->next->p;

        float dist = sqrt((visPair.p1.x - visPair.p2.x) * (visPair.p1.x - visPair.p2.x) + (visPair.p1.y - visPair.p2.y) * (visPair.p1.y - visPair.p2.y));
        visPair.dist = dist;

        visibilityGraph.push_back(visPair);
      }

      curr = curr->next;
    }
    chEdgeIntersect = chIntersect(curr->p, head->p, i);

    if (chEdgeIntersect == 0) {
      visPair.p1 = curr->p;
      visPair.p2 = head->p;

      float dist = sqrt((visPair.p1.x - visPair.p2.x) * (visPair.p1.x - visPair.p2.x) + (visPair.p1.y - visPair.p2.y) * (visPair.p1.y - visPair.p2.y));
      visPair.dist = dist;

      visibilityGraph.push_back(visPair);
    }

    /****** convex hull ends *******/

    for (int j = 0; j < obstacles.at(i).size(); j++) { // for each point of the polygon
      visPair.p1 = obstacles.at(i).at(j);

      for (int k = 0; k < obstacles.size(); k++) { // for each obstacle polygon
        if (i == k) continue; // skip if the point is on the same polygon being checked

        // get visible vertices from current point p1 for the current obstacle k
        visVertices = getVisVertices(visPair.p1, obstacles.at(k));

        // add each pair to visibility graph vector
        for (int l = 0; l < visVertices.size(); l++) {
          visPair.p2 = visVertices.at(l);

          // compute distance for the shortest path problem later
          float dist = sqrt((visPair.p1.x - visPair.p2.x) * (visPair.p1.x - visPair.p2.x) + (visPair.p1.y - visPair.p2.y) * (visPair.p1.y - visPair.p2.y));
          visPair.dist = dist;

          visibilityGraph.push_back(visPair);
        }
        visVertices.clear();
      }

      /************ compute any of the four boundary points that are visible from p1**************/
      // 1. make boundary a polygon
      vector<point2D> nonObstacle;
      point2D leftTop, rightTop, leftBottom, rightBottom;
      leftTop.x = leftBottom.x = leftBottom.y = rightBottom.y = 0;
      rightTop.x = rightBottom.x = leftTop.y = rightTop.y = WINDOWSIZE;

      nonObstacle.push_back(leftBottom);
      nonObstacle.push_back(leftTop);
      nonObstacle.push_back(rightTop);
      nonObstacle.push_back(rightBottom);

      // 2. get visible vertices of the boundary polygon
      visVertices = getVisVertices(visPair.p1, nonObstacle);

      // 3. add the pairs to the vector
      for (int l = 0; l < visVertices.size(); l++) {
        visPair.p2 = visVertices.at(l);
        
        // compute distance for the shortest path problem later
        float dist = sqrt((visPair.p1.x - visPair.p2.x) * (visPair.p1.x - visPair.p2.x) + (visPair.p1.y - visPair.p2.y) * (visPair.p1.y - visPair.p2.y));
        visPair.dist = dist;

        visibilityGraph.push_back(visPair);
      }

      /**********compute the visible points for the start and end point*********/
      nonObstacle.clear();
      visVertices.clear();

      // puth the start and end point into a vector and use getVisVertices function

      nonObstacle.push_back(start);
      nonObstacle.push_back(end);

      visVertices = getVisVertices(visPair.p1, nonObstacle);

      // add whichever point is visible
      for (int l = 0; l < visVertices.size(); l++) {
        visPair.p2 = visVertices.at(l);
        // compute distance for the shortest path problem later
        float dist = sqrt((visPair.p1.x - visPair.p2.x) * (visPair.p1.x - visPair.p2.x) + (visPair.p1.y - visPair.p2.y) * (visPair.p1.y - visPair.p2.y));
        visPair.dist = dist;

        visibilityGraph.push_back(visPair);
      }      
    }
  }
}
/****************** computing visibility graph ends here ********************/

/****************** finding shortest path ****************************/
// get index of the point in vertices array
int getIndex(point2D p) {
  for (int i = 0; i < vertices.size(); i++) {
    if (p.x == vertices.at(i).x && p.y == vertices.at(i).y) {
      return i;
    }
  }
  return -1;
}

// build adjacency lists for djikstra's to use
// in order to use indices, create the array of vertices and the array of according neighbor vertices, in the same order
void buildLists() {
  point2D a, b, c, d, start, end;
  a.x = 0;
  a.y = 0;
  b.x = 0;
  b.y = WINDOWSIZE;
  c.x = WINDOWSIZE;
  c.y = 0;
  d.x = WINDOWSIZE;
  d.y = WINDOWSIZE;
  start.x = start_x;
  start.y = start_y;
  end.x = end_x;
  end.y = end_y;

  vector<point2D> neighbors;

  // insert the vertices
  vertices.push_back(start);
  vertices.push_back(a);
  vertices.push_back(b);
  vertices.push_back(c);
  vertices.push_back(d);
  vertices.push_back(end);
  for (int i = 0; i < obstacles.size(); i++) {
    for (int j = 0; j < obstacles.at(i).size(); j++) {
      point2D vertex;
      vertex = obstacles.at(i).at(j);
      vertices.push_back(vertex);
    }
  }
  
  // go through the neighbors from each vertex, get index and weight for each neighbor vertex
  // put the vector of neighbor vertices into the adjacency list

  // neighbors of the starting point
  for (int k = 0; k < visibilityGraph.size(); k++) {
    if (visibilityGraph.at(k).p1.x == start.x && visibilityGraph.at(k).p1.y == start.y) {
      visibilityGraph.at(k).p2.index = getIndex(visibilityGraph.at(k).p2);
      visibilityGraph.at(k).p2.weight = visibilityGraph.at(k).dist;
      neighbors.push_back(visibilityGraph.at(k).p2);
    } else if (visibilityGraph.at(k).p2.x == start.x && visibilityGraph.at(k).p2.y == start.y) {
      visibilityGraph.at(k).p1.index = getIndex(visibilityGraph.at(k).p1);
      visibilityGraph.at(k).p1.weight = visibilityGraph.at(k).dist;
      neighbors.push_back(visibilityGraph.at(k).p1);
    }
  }
  adjacencyList.push_back(neighbors);
  neighbors.clear();

  for (int k = 0; k < visibilityGraph.size(); k++) {
    if (visibilityGraph.at(k).p1.x == a.x && visibilityGraph.at(k).p1.y == a.y) {
      visibilityGraph.at(k).p2.index = getIndex(visibilityGraph.at(k).p2);
      visibilityGraph.at(k).p2.weight = visibilityGraph.at(k).dist;
      neighbors.push_back(visibilityGraph.at(k).p2);
    } else if (visibilityGraph.at(k).p2.x == a.x && visibilityGraph.at(k).p2.y == a.y) {
      visibilityGraph.at(k).p1.index = getIndex(visibilityGraph.at(k).p1);
      visibilityGraph.at(k).p1.weight = visibilityGraph.at(k).dist;
      neighbors.push_back(visibilityGraph.at(k).p1);
    }
  }
  adjacencyList.push_back(neighbors);
  neighbors.clear();

  for (int k = 0; k < visibilityGraph.size(); k++) {
    if (visibilityGraph.at(k).p1.x == b.x && visibilityGraph.at(k).p1.y == b.y) {
      visibilityGraph.at(k).p2.index = getIndex(visibilityGraph.at(k).p2);
      visibilityGraph.at(k).p2.weight = visibilityGraph.at(k).dist;
      neighbors.push_back(visibilityGraph.at(k).p2);
    } else if (visibilityGraph.at(k).p2.x == b.x && visibilityGraph.at(k).p2.y == b.y) {
      visibilityGraph.at(k).p1.index = getIndex(visibilityGraph.at(k).p1);
      visibilityGraph.at(k).p1.weight = visibilityGraph.at(k).dist;
      neighbors.push_back(visibilityGraph.at(k).p1);
    }
  }
  adjacencyList.push_back(neighbors);
  neighbors.clear();

  // 
  for (int k = 0; k < visibilityGraph.size(); k++) {
    if (visibilityGraph.at(k).p1.x == c.x && visibilityGraph.at(k).p1.y == c.y) {
      visibilityGraph.at(k).p2.index = getIndex(visibilityGraph.at(k).p2);
      visibilityGraph.at(k).p2.weight = visibilityGraph.at(k).dist;
      neighbors.push_back(visibilityGraph.at(k).p2);
    } else if (visibilityGraph.at(k).p2.x == c.x && visibilityGraph.at(k).p2.y == c.y) {
      visibilityGraph.at(k).p1.index = getIndex(visibilityGraph.at(k).p1);
      visibilityGraph.at(k).p1.weight = visibilityGraph.at(k).dist;
      neighbors.push_back(visibilityGraph.at(k).p1);
    }
  }
  adjacencyList.push_back(neighbors);

  neighbors.clear();

  
  for (int k = 0; k < visibilityGraph.size(); k++) {
    if (visibilityGraph.at(k).p1.x == d.x && visibilityGraph.at(k).p1.y == d.y) {
      visibilityGraph.at(k).p2.index = getIndex(visibilityGraph.at(k).p2);
      visibilityGraph.at(k).p2.weight = visibilityGraph.at(k).dist;
      neighbors.push_back(visibilityGraph.at(k).p2);
    } else if (visibilityGraph.at(k).p2.x == d.x && visibilityGraph.at(k).p2.y == d.y) {
      visibilityGraph.at(k).p1.index = getIndex(visibilityGraph.at(k).p1);
      visibilityGraph.at(k).p1.weight = visibilityGraph.at(k).dist;
      neighbors.push_back(visibilityGraph.at(k).p1);
    }
  }
  adjacencyList.push_back(neighbors);
  neighbors.clear();

  for (int k = 0; k < visibilityGraph.size(); k++) {
    if (visibilityGraph.at(k).p1.x == end.x && visibilityGraph.at(k).p1.y == end.y) {
      visibilityGraph.at(k).p2.index = getIndex(visibilityGraph.at(k).p2);
      visibilityGraph.at(k).p2.weight = visibilityGraph.at(k).dist;
      neighbors.push_back(visibilityGraph.at(k).p2);
    } else if (visibilityGraph.at(k).p2.x == end.x && visibilityGraph.at(k).p2.y == end.y) {
      visibilityGraph.at(k).p1.index = getIndex(visibilityGraph.at(k).p1);
      visibilityGraph.at(k).p1.weight = visibilityGraph.at(k).dist;
      neighbors.push_back(visibilityGraph.at(k).p1);
    }
  }
  adjacencyList.push_back(neighbors);
  neighbors.clear();

  // go through the same order as the array of vertices
  for (int i = 0; i < obstacles.size(); i++) {
    for (int j = 0; j < obstacles.at(i).size(); j++) {
      point2D vertex;
      vertex = obstacles.at(i).at(j);

      for (int k = 0; k < visibilityGraph.size(); k++) {
        if (visibilityGraph.at(k).p1.x == vertex.x && visibilityGraph.at(k).p1.y == vertex.y) {
          visibilityGraph.at(k).p2.index = getIndex(visibilityGraph.at(k).p2);
          visibilityGraph.at(k).p2.weight = visibilityGraph.at(k).dist;
          neighbors.push_back(visibilityGraph.at(k).p2);
        } else if (visibilityGraph.at(k).p2.x == vertex.x && visibilityGraph.at(k).p2.y == vertex.y) {
          visibilityGraph.at(k).p1.index = getIndex(visibilityGraph.at(k).p1);
          visibilityGraph.at(k).p1.weight = visibilityGraph.at(k).dist;
          neighbors.push_back(visibilityGraph.at(k).p1);
        }
      }
      adjacencyList.push_back(neighbors);
      neighbors.clear();
    }
  }
}

// backtrack and return the vector of the points on the path to the given vertex
vector<point2D> getSPTo(int index, vector<int> parent) {
  vector<point2D> thisPath;
  for ( ; index != -1; index = parent[index]) {
    thisPath.insert(thisPath.begin(), vertices[index]);
  }
  return thisPath;
}

// dijkstra's algorithm using indices of each vertex and a set as a "queue" structure
void findPath(point2D s) {
  vertices.clear();
  adjacencyList.clear();
  path.clear();

  buildLists();

  // vector storing the current distance value, or status, of each vertex
  vector<float> minDistTo;
  minDistTo.clear();
  minDistTo.resize(adjacencyList.size(), 99999999.999); // initialized as a large integer
  // starting point initialized to 0
  minDistTo.at(0) = 0.0;

  // vector to save the index of previous vertex in the shortest paths to each vertex
  vector<int> parent;
  parent.clear();
  parent.resize(adjacencyList.size(), -1);

  // queue as a set of pair<weight, index>
  set<pair<float, int> > queue;
  queue.insert(make_pair(0.0, 0));

  while (!queue.empty()) {
    // queue.begin() is an element with the minimum weight due to the insert property of set
    float dist = queue.begin()->first;
    // call current vertex u
    int uIndex = queue.begin()->second;
    queue.erase(queue.begin());

    // all neighbors of the current vertex
    const vector<point2D> &uNeighbors = adjacencyList.at(uIndex);
    for (vector<point2D>::const_iterator it = uNeighbors.begin(); it != uNeighbors.end(); it++) {
      // neighbor vertex v, looking at the edge from u to v
      int vIndex = it->index;
      float weight = it->weight;
      
      // compare if the path through u to v is smaller than the current distance status of v
      // status of u + weight of the edge uv < current minimum distance to v
      float distThroughU = dist + weight;
      if (distThroughU - minDistTo[vIndex] < 0.0000001) {
        // erase the current pair and insert a new one with updated weight
        queue.erase(make_pair(minDistTo[vIndex], vIndex));

        minDistTo[vIndex] = distThroughU;
        parent[vIndex] = uIndex;
        queue.insert(make_pair(minDistTo[vIndex], vIndex));
      }
    }
  }

  // the goal is always at index 5
  path = getSPTo(5, parent);
}

/*****************finding shortest path ends here*******************/

  /* ****************************** */
void display(void) {

  glClear(GL_COLOR_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW); 
    glLoadIdentity(); //clear the matrix


    /* The default GL window is [-1,1]x[-1,1] with the origin in the
       center. 
       
       Our system of coordinates (in which we generate our points) is
       (0,0) to (WINSIZE,WINSIZE), with the origin in the lower left
       corner.
       
       We need to map the points to [-1,1] x [-1,1]  
       
       Assume we are the local coordinate system. 
       
       First we scale down to [0,2] x [0,2] */ 
      glScalef(2.0/WINDOWSIZE, 2.0/WINDOWSIZE, 1.0);  
     /* Then we translate so the local origin goes in the middle of teh
       window to (-WINDOWSIZE/2, -WINDOWSIZE/2) */
      glTranslatef(-WINDOWSIZE/2, -WINDOWSIZE/2, 0); 

    //now we draw in our local coordinate system (0,0) to
    //(WINSIZE,WINSIZE), with the origin in the lower left corner.

      draw_polygon(poly, white);
      draw_circle(start_x, start_y, cyan);
      draw_circle(end_x, end_y, yellow);

      // draw the obstacles
      draw_obstacles();

    //draw a circle where the mouse was last clicked. Note that this
    //point is stored as a global variable and is modified by the mouse
    //handler function

      draw_circle(mouse.x, mouse.y, red);

      draw_VG(red);

      draw_path(yellow);

    /* execute the drawing commands */
       glFlush();
}



  /* ****************************** */
  void keypress(unsigned char key, int x, int y) {
    switch(key) {
      case 'q':
      exit(0);
      break;

    //expected behaviour: press 's', then click on the points you
    //want, and press 'e' when you're done. the points will be saved
    //in global 'poly'

      case 's': // start entering polygons
      poly.clear();
      visibilityGraph.clear();
      path.clear();
      mouse_x = mouse_y = 0; 
      poly_init_mode = 1;
      glutPostRedisplay();
      break;

      case 'e': // end entering polygons and add the polygon to obstacles
      poly_init_mode=0;
      obstacles.push_back(poly);
      glutPostRedisplay();
      break; 

      case 'v': // compute visibility graph
      // unfinished polygons disappear
      poly_init_mode = 0;
      poly.clear();
      
      computeVG();
      glutPostRedisplay();
      break;

      case 'p': // find path using Dijkstra's algorithm
      poly_init_mode = 0;
      poly.clear();

      point2D start;
      start.x = start_x;
      start.y = start_y;

      if (visibilityGraph.size() == 0) {
        computeVG();
      }
      findPath(start);
      glutPostRedisplay();
      break;

      case 'i': // initialize the program
      poly.clear();
      poly_init_mode = 0;
      obstacles.clear();
      visibilityGraph.clear();
      vertices.clear();
      adjacencyList.clear();
      path.clear();
      // PQ_delete(pq);
      mouse_x = mouse_y = 0;
      glutPostRedisplay();
    }
  }

  void timerfunc() {

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

  int compY (const void *a, const void *b) {
    int l = ((point2D *)a)->y;
    int r = ((point2D *)b)->y;
    return (l-r);
  }

  int compX (const void *a, const void *b) {
    int l = ((point2D *)a)->x;
    int r = ((point2D *)b)->x;
    return (l-r);
  }


