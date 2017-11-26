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
  #include "pqheap.h"


  #include <stdlib.h>
  #include <stdio.h>
  #include <math.h>
  #include <assert.h>
  #include <cmath>
  #include <unistd.h>

  #ifdef __APPLE__
  #include <GLUT/glut.h>
  #else
  #include <GL/glut.h>
  #endif


  #include <vector> 
  #include <algorithm>
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
  const int WINDOWSIZE = 200;
  const int ROBOTWIDTH = 10;
  const int ROBOTHEIGHT = 20;
  const int scale = 1;

  //the current polygon 
  vector<point2D>  poly;

  //move var
  int moveInc = 0;
  int moveYet = 0;

  // obstacles as a vector of polygons
  vector<vector<point2D> > obstacles;
  // free space pre processed as a 3 dimensional vector
  // first element of the pair: 1 if free, 0 if not free
  // second element of the pair: 1 if visited, 0 if not visited
  vector<vector<vector<pair<int, int> > > > freeSpace;
  vector<vector<vector<elemType> > > parentMap;

  // vector to save the path from start to goal
  vector<elemType> path;

  point2D mouse;

  // starting point and goal for the robot
  int start_x = -1, start_y = -1, end_x = -1, end_y = -1;

/*
  int moveInc = 0;
  int xDir = 1;
  int yDir = 1;
*/

  //coordinates of last mouse click
  double mouse_x=-10, mouse_y=-10; 
  //initialized to a point outside the window

  //when this is 1, then clicking the mouse results in those points being stored in poly
  int poly_init_mode = 0; 

  // if the int is 1, clicking the mouse sets the starting point or the goal for the robot
  int start_init_mode = 0;
  int goal_init_mode = 0;

  vector<point2D> globot;

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

      if (start_init_mode == 1) {
        start_x = mouse_x;
        start_y = mouse_y;
        printf("strting point (%d, %d)\n", start_x, start_y);
      }

      if (goal_init_mode == 1) {
        end_x = mouse_x;
        end_y = mouse_y;
        printf("goal (%d, %d)\n", end_x, end_y);
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

  void initialize_freeSpace() {
    freeSpace.clear();

    freeSpace.resize(WINDOWSIZE);
    for (int x = 0; x < WINDOWSIZE; x++) {
      freeSpace.at(x).resize(WINDOWSIZE);
      for (int y = 0; y < WINDOWSIZE; y++) {
        freeSpace.at(x).at(y).resize(180/scale, make_pair(0, 0));
      }
    }
  }

  void initialize_parentMap() {
    parentMap.clear();

    elemType init;
    init.dist = init.x = init.y = init.angle = -1;
    // init.parent=NULL;

    parentMap.resize(WINDOWSIZE);
    for (int x = 0; x < WINDOWSIZE; x++) {
      parentMap.at(x).resize(WINDOWSIZE);
      for (int y = 0; y < WINDOWSIZE; y++) {
        parentMap.at(x).at(y).resize(180 / scale, init);
      }
    }
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






// return the four vertices of the robot with the center at (x, y) at a certain angle
vector<point2D> getRobot(int angle, int x, int y) {
  vector<point2D> robot;

  // four vertices before rotation
  point2D topLeft, topRight, bottomLeft, bottomRight;

  topLeft.x = x - ROBOTWIDTH / 2;
  topLeft.y = y + ROBOTHEIGHT / 2;
  topRight.x = x + ROBOTWIDTH / 2;
  topRight.y = y + ROBOTHEIGHT / 2;
  bottomLeft.x = x - ROBOTWIDTH / 2;
  bottomLeft.y = y - ROBOTHEIGHT / 2;
  bottomRight.x = x + ROBOTWIDTH / 2;
  bottomRight.y = y - ROBOTHEIGHT / 2;

  // if no rotation, return the four vertices above
  if (angle == 0) {
    robot.push_back(topLeft);
    robot.push_back(topRight);
    robot.push_back(bottomRight);
    robot.push_back(bottomLeft);
    return robot;
  }

  // vertices after rotation
  point2D v1, v2, v3, v4;

  // degree to radian
  double angleR = (double) angle * (PI / 180.0);

  // new_x = center_x + (before_x - center_x) * cos(angle) - (before_y - center_y) * sin(angle)
  // new_y = center_y + (before_x - center_x) * sin(angle) + (before_y - center_y) * cos(angle)
  v1.x = (int) (x + (topLeft.x - x) * cos(angleR) - (topLeft.y - y) * sin(angleR));
  v1.y = (int) (y + (topLeft.x - x) * sin(angleR) + (topLeft.y - y) * cos(angleR));
  v2.x = (int) (x + (topRight.x - x) * cos(angleR) - (topRight.y - y) * sin(angleR));
  v2.y = (int) (y + (topRight.x - x) * sin(angleR) + (topRight.y - y) * cos(angleR));
  v3.x = (int) (x + (bottomRight.x - x) * cos(angleR) - (bottomRight.y - y) * sin(angleR));
  v3.y = (int) (y + (bottomRight.x - x) * sin(angleR) + (bottomRight.y - y) * cos(angleR));
  v4.x = (int) (x + (bottomLeft.x - x) * cos(angleR) - (bottomLeft.y - y) * sin(angleR));
  v4.y = (int) (y + (bottomLeft.x - x) * sin(angleR) + (bottomLeft.y - y) * cos(angleR));

  robot.push_back(v1);
  robot.push_back(v2);
  robot.push_back(v3);
  robot.push_back(v4);
  return robot;
}

// if placing the robot at (angle, x, y) will cause no collision
bool isFree(int angle, int x, int y){
  // get the four vertices of the robot
  vector<point2D> robot;
  robot = getRobot(angle, x, y);


  // if the robot goes outside the boundary
  for (int i = 0; i < robot.size(); i++) {
    if (robot.at(i).x < 0 || robot.at(i).x > WINDOWSIZE ||
      robot.at(i).y < 0 || robot.at(i).y > WINDOWSIZE) {
      return 0;
    }
  }

  // check if any of the four edges intersect any obstacle
  vector<point2D> obstacle;
  
  // check for each obstacle
  for (int i = 0; i < obstacles.size(); i++) {
    obstacle = obstacles.at(i);

    // return 0 if the robot and the obstacle intersects or if one is inside the other
    if (polygonIntersect(robot, obstacle) == 1 || polygonInside(robot, obstacle) == 1 ||
      polygonInside(obstacle, robot) == 1) return 0;
  }

  return 1;
}

// go through each angle and each pixel and determine if placing the robot at that point is collision free
void getFreeSpace() {
  for (int angle = 0; angle < 180/scale; angle++) {
    for (int x = 0; x < WINDOWSIZE; x++) {
      for (int y = 0; y < WINDOWSIZE; y++) {
        if (isFree(angle * scale, x, y) == 1) {
          freeSpace.at(x).at(y).at(angle).first = 1;
        } else {
          freeSpace.at(x).at(y).at(angle).first = 0;
        }
      }
    }
  }
}


// returns the sum of euclidean distance from the start to (x, y) and from (x, y) to goal
float getDist(int x, int y) {
  return sqrt((start_x - x) * (start_x - x) + (start_y - y) * (start_y - y)) + 2 * sqrt((end_x - x) * (end_x - x) + (end_y - y) * (end_y - y));
}

vector<elemType> getSuccessors(elemType *s) {
  int newAngle, newX, newY;
  vector<elemType> successors;

  //printf("successors for (%d, %d, %d)\n", s->x, s->y, s->angle);
  for (int i = -1; i <= 1; i++) {
    // skip if out of boundary
    if (s->x + i < 0 || s->x + i >= WINDOWSIZE) continue;
    newX = s->x + i;
    
    for (int j = -1; j <= 1; j++) {
      // skip if out of boundary
      if (s->y + j < 0 || s->y + j >= WINDOWSIZE) continue;
      newY = s->y + j;
    
      for (int k = -1; k <= 1; k++) {
        // skip (+0, +0, +?)
        if (i == 0 && j == 0 && k == 0) continue;
        // skip if out of boundary
        if (s->angle + k < 0 || s->angle + k >= 180) continue;
        newAngle = s->angle +k;

        // add the new state to the vector of successors
        elemType newState;
        newState.angle = newAngle;
        newState.x = newX;
        newState.y = newY;
        newState.dist = getDist(newX, newY);
        
        successors.push_back(newState);
      }
    }
  }

  return successors;
}

void renderPath(elemType s) {
  /*
  elemType *curr = s;
  printf("test\n");
  while (curr->parent != NULL) {
    thisPath.insert(thisPath.begin(), *curr);
    printf("before (%f, %d, %d, %d)\n", curr->parent->dist, curr->parent->x, curr->parent->y, curr->parent->angle);
    curr = curr->parent;
    printf("    after (%f, %d, %d, %d)\n", curr->dist, curr->x, curr->y, curr->angle);
  }
  */

  for (int angle = 0; angle < 180/scale; angle++) {
    for (int x = 0; x < WINDOWSIZE; x++) {
      for (int y = 0; y < WINDOWSIZE; y++) {
        freeSpace.at(x).at(y).at(angle).second = 0;
      }
    }
  }


  elemType curr = s;
  while (curr.x != start_x && curr.y != start_y) {
    path.insert(path.begin(), curr);
    curr = parentMap.at(curr.x).at(curr.y).at(curr.angle/scale);
  }
  elemType start;
  start.x = start_x;
  start.y = start_y;
  start.angle = curr.angle;
  path.insert(path.begin(), start);
  return;
}


void findPath() {
  if (start_x < 0 || start_y < 0) {
    printf("Press R to set the starting point for the robot\n");
    return;
  }

  if (end_x < 0 || end_y < 0) {
    printf("Press G to set the goal for the robot\n");
    return;
  }

  if (obstacles.size() == 0) {
    printf("Press S to start making obstacles\n");
    return;
  }

  if (!isObstaclesSimple()) {
    printf("The obstacles are not simple polygons. Press i to initialize the scene and s to draw new obstacles.\n");
    return;
  }

  initialize_parentMap();
  moveInc = moveYet = 0;

  PQueue *pq;
  pq = PQ_initialize();

  elemType start;
  start.angle = 0;
  start.x = start_x;
  start.y = start_y;
  start.dist = getDist(start.x, start.y);
  // start.parent = NULL;
  PQ_insert(pq, start);

  // if (start.parent == NULL) printf("start parent is null\n");
  
  elemType min, successor, goal;
  vector<elemType> successors;
  int reachedGoal = 0;
  while (reachedGoal != 1) {
    if (PQ_isEmpty(pq)) {
      printf("Couldn't find the path. Press i to initialize the program\n");
      // do something to delete the current path
      break;
    }
    
    assert(PQ_extractMin(pq, &min));
    successors = getSuccessors(&min);
/*
    printf("extracted (%f, %d, %d, %d)\n", min.dist, min.x, min.y, min.angle);
    for (int i = 0; i < successors.size() - 1; i++) {
      printf("successor (%f, %d, %d, %d) has parent (%f, %d, %d, %d)\n", successors.at(i).dist, successors.at(i).x, successors.at(i).y, successors.at(i).angle, (successors.at(i).parent)->dist, (successors.at(i).parent)->x, (successors.at(i).parent)->y, (successors.at(i).parent)->angle);
    }
*/
/*
    // insert the successors that are free without rotation first
    int rotationFreeCount = 0;
    for (vector<elemType>::iterator it = successors.begin(); it != successors.end(); ++it) {
      if (it->x == end_x && it->y == end_y) {
        reachedGoal = 1;
        goal = *it;
        printf("found goal!\n");
        rotationFreeCount = -1;
        break;
      }

      successor = *it;

      // if the program has already looked at the successor state, skip
      if (freeSpace.at(successor.x).at(successor.y).at(successor.angle).second == 1) continue;

      freeSpace.at(successor.x).at(successor.y).at(successor.angle).second = 1;

      // if the successor state (x, y, angle) is free
      if (freeSpace.at(successor.x).at(successor.y).at(successor.angle).first == 1) {
        // if the successor state is free without rotation
        if (successor.angle == min.angle) {
          rotationFreeCount++;
          PQ_insert(pq, successor);
        }
      }
    }
    // skip to next state in the queue if any free successor without rotation is found
    if (rotationFreeCount != 0) continue;
*/
    // add the rest of the free successors to the queue
    for (vector<elemType>::iterator i = successors.begin(); i != successors.end(); ++i) {
      if (i->x == end_x && i->y == end_y) {
        reachedGoal = 1;
        goal = *i;
        parentMap.at(goal.x).at(goal.y).at(goal.angle / scale) = min;
        //goal.parent = &min;
        //printf("goal parent (%f, %d, %d, %d), whose parent is (%f, %d, %d, %d)\n", goal.parent->dist, goal.parent->x, goal.parent->y, goal.parent->angle, goal.parent->parent->dist, goal.parent->parent->x, goal.parent->parent->y, goal.parent->parent->angle);
        break;
      }

      successor = *i;

      
      //printf("parent (%d, %d, %d)\n", (successor.parent)->x, (successor.parent)->y, (successor.parent)->angle);
      // if the program has already looked at the successor state, skip
      if (freeSpace.at(successor.x).at(successor.y).at(successor.angle / scale).second == 1) continue;
      // mark as visited
      freeSpace.at(successor.x).at(successor.y).at(successor.angle / scale).second = 1;

      parentMap.at(successor.x).at(successor.y).at(successor.angle / scale) = min;      

      // if the successor state (angle, x, y) is free, add it to the queue
      if (freeSpace.at(successor.x).at(successor.y).at(successor.angle / scale).first == 1) {
        // successor.parent = &min;
        //printf("setting parent (%f, %d, %d, %d)\n", min.dist, min.x, min.y, min.angle);
        PQ_insert(pq, successor);
      }
    }
  }

  PQ_delete(pq);

  if (goal.x == end_x && goal.y == end_y) {
    renderPath(goal); // have to use *(whatever.parent) to get the parent elemType
    moveYet = 1;
  }
}


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

      if (moveYet == 1) {
        draw_polygon(globot, green);
      }

    //draw a circle where the mouse was last clicked. Note that this
    //point is stored as a global variable and is modified by the mouse
    //handler function

      draw_circle(mouse.x, mouse.y, red);

      //draw_path(yellow);
      
      

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
      case 'r':
      poly.clear();
      poly_init_mode = 0;
      goal_init_mode = 0;
      mouse_x = mouse_y = -1;
      start_init_mode = 1;

      break;

      case 'g':
      poly.clear();
      poly_init_mode = 0;
      start_init_mode = 0;
      mouse_x = mouse_y = -1;
      goal_init_mode = 1;
      break;

      case 's': // start entering polygons
      start_init_mode = 0;
      goal_init_mode = 0;
      poly.clear();
      path.clear();
      mouse_x = mouse_y = 0;
      poly_init_mode = 1;
      glutPostRedisplay();
      break;

      case 'e': // end entering polygons and add the polygon to obstacles
      poly_init_mode=0;
      start_init_mode = 0;
      goal_init_mode = 0;
      if (poly.size() > 2) {
        obstacles.push_back(poly);
      }
      glutPostRedisplay();
      break; 

/*
      case 'v': // compute visibility graph
      // unfinished polygons disappear
      poly_init_mode = 0;
      start_init_mode = 0;
      poly.clear();
      
      computeVG();
      glutPostRedisplay();
      break;

      */

      case 'f':
      poly_init_mode=0;
      start_init_mode = 0;
      goal_init_mode = 0;
      initialize_freeSpace();
      getFreeSpace();
      break;

      case 'p': // find path using Dijkstra's algorithm
      poly_init_mode = 0;
      start_init_mode = 0;
      goal_init_mode = 0;
      poly.clear();
      path.clear();
      globot.clear();
      moveYet = 0;
      findPath();
      glutPostRedisplay();
      break;

      case 'i': // initialize the program
      poly.clear();
      poly_init_mode = 0;
      start_init_mode = 0;
      goal_init_mode = 0;
      obstacles.clear();
      path.clear();
      globot.clear();
      moveYet = 0;
      mouse_x = mouse_y = -1;
      start_x = start_y = -1;
      end_x = end_y = -1;
      glutPostRedisplay();
      break;
    }
  }

  void timerfunc() {

    if (moveYet == 1) {
      if (moveInc < path.size()) {
        //vector<point2D> robot;
        printf("Path %d (angle = %d, x = %d, y = %d)\n", moveInc+1, path.at(moveInc).angle, path.at(moveInc).x, path.at(moveInc).y);
        globot = getRobot(path.at(moveInc).angle, path.at(moveInc).x, path.at(moveInc).y);
        
      }
      moveInc += 1;;
    }
    glutPostRedisplay();
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


