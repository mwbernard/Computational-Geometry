  /* mouse2.c 

  Jeonguk Choi and Martin Bernard

  The program calculates and displays visibility polygon of the point named "guard."
  When s is pressed, mouse clicks create a polygon.
  When e is pressed, adding polygon stops.
  When c is pressed, computes and displays the visibility polygon of the last mouse click point.
  When m is pressed, continuously computes and displays the visibility polygon of the moving point.
  */

  #include "geom.h"


  #include <stdlib.h>
  #include <stdio.h>
  #include <math.h>
  #include <assert.h>
  #include <cmath>

  #ifdef __APPLE__
  #include <GLUT/glut.h>
  #else
  #include <GL/glut.h>
  #endif


  #include <vector> 

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

  void initialize_polygon(); 
  void print_polygon(vector<point2D> poly); 
  point2D getSecondMin(point2D p, vector<point2D> points);
  vector<point2D> findIntersections(point2D p, point2D vertex);
  void draw_triangles(vector<point2D> points, point2D p);
  bool isSimple();
  void visPolygon();
  point2D intersect(point2D a, point2D b, point2D c, point2D d);
  char isInside(point2D p);
  vector<point2D> getVisPoints(vector<point2D> poly, point2D p);
  int compY (const void *a, const void *b);
  int compX (const void *a, const void *b);

  int left (point2D a, point2D b, point2D c);
  int signed_area2D(point2D a, point2D b, point2D c);
  int collinear(point2D p, point2D q, point2D r);
  int collinear2(point2D p, point2D q, point2D r);
  double distance(point2D a, point2D b);

  /* our coordinate system is (0,0) x (WINDOWSIZE,WINDOWSIZE) where the
     origin is the lower left corner */


  /* global variables */
  const int WINDOWSIZE = 750; 

  //the current polygon 
  vector<point2D>  poly;
  vector<point2D> visPoints;
  point2D guard;
  int moveBool = 0;
  int xDir = 1;
  int yDir = 1;

  //coordinates of last mouse click
  double mouse_x=-10, mouse_y=-10; 
  //initialized to a point outside the window

  //when this is 1, then clicking the mouse results in those points being stored in poly
  int poly_init_mode = 0; 



  void orderPoly(){

    point2D p1, p2;
    int sum = 0;

    for (int i = 0; i < poly.size() - 1; i++) {
      p1 = poly.at(i);
      p2 = poly.at(i+1);
      sum += (p2.x - p1.x)*(p2.y + p1.y);
    }

    p1 = poly.at(poly.size()-1);
    p2 = poly.at(0);
    sum += (p2.x - p1.x)*(p2.y + p1.y);

    if (sum < 0){ 
      reverse(poly.begin(), poly.end());
    }
  }

  void draw_circle(double x, double y){
    
    glColor3fv(red);   
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
      visPoints.clear();

      mouse_x = x;
      mouse_y = y;
      //(x,y) are in wndow coordinates, where the origin is in the upper
      //left corner; our reference system has the origin in lower left
      //corner, this means we have to reflect y
      mouse_y = WINDOWSIZE - mouse_y; 
      if (poly_init_mode ==1) {
        point2D p = {mouse_x, mouse_y, -1}; 
        poly.push_back(p);
      }
    }
    
    guard.x = (int) mouse_x;
    guard.y = (int) mouse_y;
    glutPostRedisplay();
  }




  /* ****************************** */
  /* initialize  polygon stored in global variable poly  */
  void initialize_polygon() {

    //clear the vector, in case something was there 
    poly.clear(); 

    int n = 10; //size of polygon 
    double rad = 100; 
    double step = 2 * M_PI / n;
    point2D p; 
    for (int i=0; i<n; i++) {

      p.x = WINDOWSIZE/2 + rad * cos (i * step); 
      p.y = WINDOWSIZE/2 + rad * sin (i * step); 


      //insert the segment in the array of segments 
      poly.push_back (p); 
    } //for i
    
  }



  /* ************************************************** */
  void print_polygon(vector<point2D> poly) {

    for (unsigned int i=0; i<poly.size()-1; i++) {
      printf("edge %d: [(%d,%d), (%d,%d)]\n",
        i, poly[i].x, poly[i].y, poly[i+1].x, poly[i+1].y);
    }
    //print last edge from last point to first point 
    int last = poly.size()-1; 
    printf("edge %d: [(%d,%d), (%d,%d)]\n",
      last, poly[last].x, poly[last].y, poly[0].x, poly[0].y);

  }
/*
  decides directions based on the velocity vectors
  doesn't really work well as of now....

  void moveGuard() {

    if (moveBool == 1) {
      
      point2D v1;
      point2D v2;

      point2D temp;
      temp.x = guard.x + 5 * xDir;
      temp.y = guard.y + 5 * yDir;

      int biggerX, smallerX, biggerY, smallerY;

      for (int i = 0; i < poly.size() - 1; i++){

        v1 = poly.at(i);
        v2 = poly.at(i + 1);

        if (v1.x > v2.x) {
          biggerX = v1.x;
          smallerX = v2.x;
        } else {
          biggerX = v2.x;
          smallerX = v1.x;
        }
        if (v1.y > v2.y) {
          biggerY = v1.y;
          smallerY = v2.y;
        } else {
          biggerY = v2.y;
          smallerY = v1.y;
        }

        // collision
        if (left(v1, v2, temp) != left(v1, v2, guard) && guard.x <= biggerX && guard.x >= smallerX && guard.y <= biggerY && guard.y >= smallerY) {
          
          double dx = (double) (v2.x - v1.x) / distance(v2, v1);
          double dy = (double) (v2.y - v1.y) / distance(v2, v1);

          // get two normal vectors for each edge
          // two sets of vectors are necessary in order to calculate the direction of the bouncing vector
          double uX1, uY1, wX1, wY1, uX2, uY2, wX2, wY2, xDir1, yDir1, xDir2, yDir2;

          // normal vector <-dy, dx>
          uX1 = ((xDir * (-1) * dy) + (yDir * dx)) * (-1) * dy;
          uY1 = ((xDir * (-1) * dy) + (yDir * dx)) * dx;

          wX1 = xDir - uX1;
          wY1 = yDir - uY1;

          xDir1 = wX1 - uX1;
          yDir1 = wY1 - uY1;

          // normal vector <dy, -dx>
          uX2 = ((xDir * dy) + (yDir * (-1) * dx)) * dy;
          uY2 = ((xDir * dy) + (yDir * (-1) * dx)) * (-1) * dx;

          wX2 = xDir - uX2;
          wY2 = yDir - uY2;

          xDir2 = wX2 - uX2;
          yDir2 = wY2 - uY2;

          point2D nextPoint1, nextPoint2;
          nextPoint1.x = guard.x + 5 * xDir1;
          nextPoint1.y = guard.y + 5 * yDir1;
          nextPoint2.x = guard.x + 5 * xDir2;
          nextPoint2.y = guard.y + 5 * yDir2;

          // decide the new x and y directions based on which of the two directions gives a point inside
          if (isInside(nextPoint1) == 'i') {
            xDir = wX1 - uX1;
            yDir = wY1 - uY1;
          } else if (isInside(nextPoint2) == 'i') {
            xDir = wX2 - uX2;
            yDir = wY2 - uY2;
          } else {
            xDir = (-1) * xDir;
            yDir = (-1) * yDir;
          }
        } 
      }

      int end = poly.size();
      v1 = poly.at(end - 1);
      v2 = poly.at(0);

      if (v1.x > v2.x) {
          biggerX = v1.x;
          smallerX = v2.x;
        } else {
          biggerX = v2.x;
          smallerX = v1.x;
        }
        if (v1.y > v2.y) {
          biggerY = v1.y;
          smallerY = v2.y;
        } else {
          biggerY = v2.y;
          smallerY = v1.y;
        }

        // collision
      if (left(v1, v2, temp) != left(v1, v2, guard) && guard.x <= biggerX && guard.x >= smallerX && guard.y <= biggerY && guard.y >= smallerY) {
        double dx = (double) (v2.x - v1.x) / distance(v2, v1);
        double dy = (double) (v2.y - v1.y) / distance(v2, v1);

        // there are two normal vectors for each edge
        // two sets of vectors are necessary in order to calculate the direction of the bouncing vector
        double uX1, uY1, wX1, wY1, uX2, uY2, wX2, wY2, xDir1, yDir1, xDir2, yDir2;

        // normal vector <-dy, dx>
        uX1 = ((xDir * (-1) * dy) + (yDir * dx)) * (-1) * dy;
        uY1 = ((xDir * (-1) * dy) + (yDir * dx)) * dx;

        wX1 = xDir - uX1;
        wY1 = yDir - uY1;

        xDir1 = wX1 - uX1;
        yDir1 = wY1 - uY1;

        // normal vector <dy, -dx>
        uX2 = ((xDir * dy) + (yDir * (-1) * dx)) * dy;
        uY2 = ((xDir * dy) + (yDir * (-1) * dx)) * (-1) * dx;

        wX2 = xDir - uX2;
        wY2 = yDir - uY2;

        xDir2 = wX2 - uX2;
        yDir2 = wY2 - uY2;

        point2D nextPoint1, nextPoint2;
        nextPoint1.x = guard.x + 5 * xDir1;
        nextPoint1.y = guard.y + 5 * yDir1;
        nextPoint2.x = guard.x + 5 * xDir2;
        nextPoint2.y = guard.y + 5 * yDir2;

        // decide the new x and y directions based on which of the two directions gives a point inside
        if (isInside(nextPoint1) == 'i') {
          xDir = wX1 - uX1;
          yDir = wY1 - uY1;
        } else if (isInside(nextPoint2) == 'i') {
          xDir = wX2 - uX2;
          yDir = wY2 - uY2;
        } else {
          xDir = (-1) * xDir;
          yDir = (-1) * yDir;
        }
      }
      
      guard.x += 2*xDir;
      guard.y += 2*yDir;
    }
  }
*/

  // randomly selects x- and y-direction of the moving point, choosing from -1, 0, and 1
  void moveGuard() {
    if (moveBool == 1) {
      // vertices of the edge
      point2D v1;
      point2D v2;

      // two steps ahead in the current direction
      // in order to check if the point is on the other side of the edge
      point2D temp;
      temp.x = guard.x + 2 * xDir;
      temp.y = guard.y + 2 * yDir;

      int biggerX, smallerX, biggerY, smallerY;

      for (int i = 0; i < poly.size(); i++){
        // for each edge
        if (i == poly.size() - 1) {
          v1 = poly.at(i);
          v2 = poly.at(0);
        } else {
          v1 = poly.at(i);
          v2 = poly.at(i + 1);
        }

        // boundary coordinates to determine if the point is within the edge
        if (v1.x > v2.x) {
          biggerX = v1.x;
          smallerX = v2.x;
        } else {
          biggerX = v2.x;
          smallerX = v1.x;
        }
        if (v1.y > v2.y) {
          biggerY = v1.y;
          smallerY = v2.y;
        } else {
          biggerY = v2.y;
          smallerY = v1.y;
        }

        // collision when the guard and the point a few steps ahead in current direction 
        // are on different sides of the current edge and within the line
        if (left(v1, v2, temp) != left(v1, v2, guard) && guard.x <= biggerX 
          && guard.x >= smallerX && guard.y <= biggerY && guard.y >= smallerY) {
          int newXDir, newYDir;
          point2D nextPoint;

          do {
            // pick a random value between -1 and 1
            newXDir = rand() % 3 - 1;
            newYDir = rand() % 3 - 1;
            
            // if both x and y directions are 0, pick a new random pair of values
            while (newXDir == 0 && newYDir == 0) {
              newXDir = rand() % 3 - 1;
              newYDir = rand() % 3 - 1;
            }

            // two steps ahead in the direction to check if it's still inside
            nextPoint.x = guard.x + 2 * newXDir;
            nextPoint.y = guard.y + 2 * newYDir;
          } while (isInside(nextPoint) != 'i'); // if the point will go outside, pick a new direction

          // set the direction to the newly assembled direction
          xDir = newXDir;
          yDir = newYDir;
        }
      }
      // move the guard's position
      guard.x = guard.x + 2 * xDir;
      guard.y = guard.y + 2 * yDir;
    }
  }


  /* ****************************** */
  int main(int argc, char** argv) {

    initialize_polygon();
    print_polygon(poly);


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
  void draw_polygon(vector<point2D> poly, GLfloat *color){

    if (poly.size() == 0) return; 

    //set color 
    glColor3fv(color);

    glLineWidth(3);  
    
    int i;
    for (i=0; i<poly.size()-1; i++) {
      glBegin(GL_LINES);
      glVertex2f(poly[i].x, poly[i].y); 
      glVertex2f(poly[i+1].x, poly[i+1].y);
      glEnd();
    }
    //render last segment between last point and forst point 
    int last=poly.size()-1; 
    glBegin(GL_LINES);
    glVertex2f(poly[last].x, poly[last].y); 
    glVertex2f(poly[0].x, poly[0].y);
    glEnd();
  }

  // insert the current vertex (because there are cases where the vertex is not found due to approximation errors)
  // remove the duplicate points in the vector
  vector<point2D> dupRemoval(point2D vertex, vector<point2D> points) { 
    
    points.push_back(vertex);

    if (points.size() < 2) return points;

    point2D p1, p2;
    for (vector<point2D>::iterator it = points.begin(); it != points.end(); ++it) {
      p1 = *it;
      for (vector<point2D>::iterator it2 = points.begin(); it2 != points.end(); ++it2) {
        if (it == it2) continue;
        p2 = *it2;

        if (p1.x - p2.x < 5 && p1.x - p2.x > -5 && p1.y - p2.y < 5 && p1.y - p2.y > -5) {
          points.erase(it2);
          --it2;
        }
      }
    }

    return points;
  }

  // get visible vertices from the point p
  // if the closest intersection between p-vertex and each edge is that vertex, it is visible
  vector<point2D> getVisVertices(point2D p) {
    vector<point2D> visVertices;

    // check each line between p and a vertex
    for(vector<point2D>::iterator it = poly.begin(); it != poly.end(); ++it) {
      point2D curr = *it;
      vector<point2D> edgeIntersections;

      // get all the intersections between p-curr and every edge
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

      if (minDistIntersection.vertex == 1 && abs(minDistIntersection.x - curr.x) < 5) {
          visVertices.push_back(minDistIntersection);
      }
    }

    return visVertices;
  }

  // calculate the second closest intersection point between each edge and the ray from p in the direction of each vertex
  point2D getExtLineVisPoint(point2D vertex, point2D p, int xExtreme, int yExtreme) {
    point2D extPoint; // extreme endpoint of the ray
    vector<point2D> extLineIntersections; 
    point2D extLineVisible; // visible point i.e. the second closest intersection
    double x_intercept, y_intercept;

    y_intercept = (double)((xExtreme - vertex.x)*(vertex.y - p.y)) / (double) (vertex.x - p.x) + vertex.y;
    x_intercept = (double)((yExtreme - vertex.y)*(vertex.x - p.x)) / (double) (vertex.y - p.y) + vertex.x;
  
    if ((x_intercept < xExtreme && xExtreme == WINDOWSIZE) || (x_intercept > xExtreme && xExtreme == 0)) {
      extPoint.x = x_intercept;
      extPoint.y = yExtreme;
      extLineIntersections = findIntersections(p, extPoint);
      extLineIntersections = dupRemoval(vertex, extLineIntersections);
      extLineVisible  = getSecondMin(p, extLineIntersections);
    }
    else if ((y_intercept < yExtreme && yExtreme == WINDOWSIZE) || (y_intercept > yExtreme && yExtreme == 0)) {
      extPoint.x = xExtreme;
      extPoint.y = y_intercept;
      extLineIntersections = findIntersections(p, extPoint);
      extLineIntersections = dupRemoval(vertex, extLineIntersections);
      extLineVisible  = getSecondMin(p, extLineIntersections);
    }
    else {
      extPoint.x = xExtreme;
      extPoint.y = yExtreme;
      extLineIntersections = findIntersections(p, extPoint);
      extLineIntersections = dupRemoval(vertex, extLineIntersections);
      extLineVisible = getSecondMin(p, extLineIntersections);
    }

    return extLineVisible;
  }

  // compute all the visible points on the polygon from p
  void getAllVisPoints(point2D p, vector<point2D> visVertices) {
    point2D extLineVisPoint; // visible point calculated from the ray from p in the direction of each vertex

    for (int i = 0; i < visVertices.size(); i++) {
      point2D prev, vertex, next;
      for (int j = 0; j < poly.size(); j++) {
        if (visVertices.at(i).x == poly.at(j).x && visVertices.at(i).y == poly.at(j).y) {
          vertex = visVertices.at(i);
          if (j == 0) {
            prev = poly.at(poly.size() - 1);
            next = poly.at(1);
          } else if (j == poly.size() - 1) {
            prev = poly.at(j - 1);
            next = poly.at(0);
          } else {
            prev = poly.at(j - 1);
            next = poly.at(j + 1);
          }
        }
      }

      // if the two outgoing edges from a visible vertex is on the same side of the line p-vertex
      if ((left(p, vertex, prev) == left(p, vertex, next)) || (left(p, vertex, prev) && collinear(p, vertex, next))
       || (left(p, vertex, next) && collinear(p, vertex, prev))) {
        vector<point2D> extLineIntersections;
        point2D extPoint;

        // different endpoint of the ray for each direction
        if ((vertex.x - p.x) > 0) { // x+
          if ((vertex.y - p.y) > 0) { // x+ y+
            extLineVisPoint = getExtLineVisPoint(vertex, p, WINDOWSIZE, WINDOWSIZE);
          }
          else if ((vertex.y - p.y) < 0) { // x+ y-
            extLineVisPoint = getExtLineVisPoint(vertex, p, WINDOWSIZE, 0);
          }
          else { // x+ y0
            extPoint.x = WINDOWSIZE;
            extPoint.y = p.y;
            extLineIntersections = findIntersections(p, extPoint);
            extLineIntersections = dupRemoval(vertex, extLineIntersections);
            extLineVisPoint  = getSecondMin(p, extLineIntersections);
          }
        } else if ((vertex.x - p.x) < 0) { // x-
          if ((vertex.y - p.y) > 0) { // x- y+
            extLineVisPoint = getExtLineVisPoint(vertex, p, 0, WINDOWSIZE);
          }
          else if ((vertex.y - p.y) < 0) { // x- y-
            extLineVisPoint = getExtLineVisPoint(vertex, p, 0, 0);  
          }
          else { // x- y0
            extPoint.x = 0;
            extPoint.y = p.y;
            extLineIntersections = findIntersections(p, extPoint);
            extLineIntersections = dupRemoval(vertex, extLineIntersections);
            extLineVisPoint  = getSecondMin(p, extLineIntersections);
          }
        } else { // x0
          if ((vertex.y - p.y) > 0) { // x0 y+
            extPoint.y = WINDOWSIZE;
            extPoint.x = p.x;
            extLineIntersections = findIntersections(p, extPoint);
            extLineIntersections = dupRemoval(vertex, extLineIntersections);
            extLineVisPoint  = getSecondMin(p, extLineIntersections);
          }
          else if ((vertex.y - p.y) < 0) { // x0 y-
            extPoint.y = 0;
            extPoint.x = p.x;
            extLineIntersections = findIntersections(p, extPoint);
            extLineIntersections = dupRemoval(vertex, extLineIntersections);
            extLineVisPoint  = getSecondMin(p, extLineIntersections);
          }
        }

        // determine the order of insertion for the new points
        // if outgoing edges are on the left, vertex first; if they are on the right, the new point first
        if (left(p, vertex, prev)) {
          visPoints.push_back(vertex);
          visPoints.push_back(extLineVisPoint);
        } else {
          visPoints.push_back(extLineVisPoint);
          visPoints.push_back(vertex);
        }
      }
      // two outgoing edges are on the different sides
      else {
        visPoints.push_back(vertex);
      }
    }
  }

  // calculates the visible polygon of a given point
  void visPolygon(){
    
    visPoints.clear();
    orderPoly();    

    // check if the polygon is simple
    if(!isSimple()) {
      printf("The polygon is not simple.\n");
      // do something to make the user draw a new polygon
      printf("Press s and draw a new polygon.\n");
    }
    
    // check if the point is inside the polygon
    if(isInside(guard) == 'o'){
      printf("The point is not inside the polygon\n");
    }

    if (isInside(guard) != 'o' && isSimple()) {
      // vector of visible vertices from p
      vector<point2D> visVertices;
      visVertices = getVisVertices(guard);

      // get all the visible points on the polygon and add them to the vector visPoints
      getAllVisPoints(guard, visVertices);

    	// draw the triangles
      draw_triangles(visPoints, guard);
      glFlush();
    }
  }

  // draws the triangles connecting p and two consecutive visible points
void draw_triangles(vector<point2D> points, point2D p) {
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  for (int i = 0; i < points.size(); i++) {
    if (i == points.size() - 1) {
      glColor3fv(cyan);
      glBegin(GL_TRIANGLES);
      glVertex2f(points.at(i).x, points.at(i).y);
      glVertex2f(points.at(0).x, points.at(0).y);
      glVertex2f(p.x, p.y);
      glEnd();
    } else {
      glColor3fv(cyan);
      glBegin(GL_TRIANGLES);
      glVertex2f(points.at(i).x, points.at(i).y);
      glVertex2f(points.at(i + 1).x, points.at(i + 1).y);
      glVertex2f(p.x, p.y);
      glEnd();
    }
  }
}

// return the second closest point from p
point2D getSecondMin(point2D p, vector<point2D> points) {
  int minDist = 99999999;
  int secondMinDist = 99999999;
  point2D min = p;
  point2D secondMin = p;

  if (points.size() == 0) {
    return p;
  } else if (points.size() == 1) {
    return points.at(0);
  }

  for (int i = 0; i < points.size(); i++) {
    int dist;
    dist = (p.x - points.at(i).x) * (p.x - points.at(i).x) +
    (p.y - points.at(i).y) * (p.y - points.at(i).y);

    if (dist < minDist) {
      secondMinDist = minDist;
      secondMin = min;
      minDist = dist;
      min = points.at(i);
    } else if (dist < secondMinDist) {
      secondMinDist = dist;
      secondMin = points.at(i);
    }
  }

  return secondMin;
}

// returns the intersections between p-vertex and each edge
vector<point2D> findIntersections(point2D p, point2D vertex) {
  vector<point2D> intersections;
  for (int i = 0; i < poly.size(); i++) {
    point2D intersection;
    point2D endpoint1, endpoint2;
    if (i == poly.size() - 1) {
      endpoint1 = poly.at(i);
      endpoint2 = poly.at(0);
      intersection = intersect(p, vertex, poly.at(i), poly.at(0));
    } else {
      endpoint1 = poly.at(i);
      endpoint2 = poly.at(i + 1);
      intersection = intersect(p, vertex, poly.at(i), poly.at(i + 1));
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
  return intersections;
}

  // Because drawing function connects the last point to the first point, there cannot be any hole
  // Therefore, the polygon is not simple only if it has any intersection in the middle of the line
bool isSimple() {
  point2D a, b, c, d;

  // go through each pair of edges and check if there is any intersection

  for (vector<point2D>::iterator it1 = poly.begin(); it1 != poly.end(); ++it1) {
    a = *it1;
    if (it1 == poly.end() - 1) {
      b = poly.at(0);
    } else {
      b = *(it1 + 1);
    }

    for (vector<point2D>::iterator it2 = poly.begin(); it2 != poly.end(); ++it2) {
      if (it1 == it2) continue;

      c = *it2;
      if (it2 == poly.end() - 1) {
        d = poly.at(0);
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
char isInside(point2D p) { 
  int i, i1;
  double x;
  int rCross = 0;
  int lCross = 0;
  bool rStrad, lStrad;
  vector<point2D> shiftedPolygon;

  for (i = 0; i < poly.size(); i++) {
    shiftedPolygon.push_back(poly.at(i));
    shiftedPolygon.at(i).x = shiftedPolygon.at(i).x - p.x;
    shiftedPolygon.at(i).y = shiftedPolygon.at(i).y - p.y;
  }

  for (i = 0; i < poly.size(); i++) {
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

    //draw a circle where the mouse was last clicked. Note that this
    //point is stored as a global variable and is modified by the mouse
    //handler function

      draw_circle(guard.x, guard.y);


      draw_triangles(visPoints, guard);
      draw_circle(guard.x, guard.y);
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

      case 's': 
      poly.clear();
      visPoints.clear();
      mouse_x = mouse_y=0; 
      poly_init_mode = 1;
      moveBool = 0;
      glutPostRedisplay();
      break; 

      case 'e': 
      poly_init_mode=0;
      moveBool = 0;
      visPoints.clear();
      guard.x = mouse_x;
      guard.y = mouse_y;
      glutPostRedisplay();
      break; 

      case 'c':
      moveBool = 0;
      visPolygon();
      glutPostRedisplay();
      break;

      case 'm':
      point2D p;
      p.x = mouse_x;
      p.y = mouse_y;
      if (poly.size() != 0 && isInside(p) == 'i') {
        moveBool = 1;
      } else {
        printf("Point is not inside the polygon");
        moveBool = 0;
      }
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

  void timerfunc() {
    if (moveBool == 1) {
      moveGuard();
      if (isInside(guard) != 'i') {
        xDir = (-1) * xDir;
        yDir = (-1) * yDir;
        guard.x = guard.x + xDir;
        guard.y = guard.y + yDir;
      }
      visPolygon();
    }
    glutPostRedisplay(); 
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

  double distance(point2D a, point2D b) {
    return sqrt(((a.x - b.x) * (a.x - b.x)) + ((a.y - b.y) * (a.y - b.y))); 
  }