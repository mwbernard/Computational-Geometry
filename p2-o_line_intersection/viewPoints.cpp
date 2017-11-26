/* view.c 

Laura Toma

What it does:  

Draws a set of horizontal and vertical line segments in the default 2D
projection. Then computes their intersections using the line sweep
algorithm, and  simulates the algorithm as it runs.

*/

#include "geom.h"
#include "rtimer.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


#include <vector>
#include <map>
#include <iostream>

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
void timerfunc(); 

void initialize_segments_random(); 
void initialize_segments_horizontal();
void initialize_segments_vertical();
void initialize_segments_overlap();
void initialize_segments_grid();
void print_segments(); 

//draws a segment s, with a desired width
void draw_segment(segment2D s, GLfloat width);
//renders the sweep line 
void draw_sweep_line();
//renders the active structure
void draw_active_structure();
//renders the intersection points 
void draw_intersection_points();


void intersection();
void getEvents();
bool comp(const sweepEvent &a, const sweepEvent &b);
vector<segment2D> rangeSearch(int y1, int y2);
void addPoints(vector<segment2D> s, int x);

/* global variables */
const int WINDOWSIZE = 500; 

int init_case = 0; 
const int NB_TEST_CASES = 5;

//NOTE: all the structures below need to be global so that they can be rendered

//current position of sweep line 
int sweep_line_x = 0; 


//number of segments requested by user 
int n; 

//the array of  segments
vector<segment2D> segments;

//the intersections points of the segments 
vector<point2D> intpoints; 

//the active structure that stores the segments intersecting the sweep line 
multimap<int, segment2D> as;

//the events 
vector<sweepEvent> events; 




void intersection() {
    as.clear();
    intpoints.clear();
    events.clear();
  getEvents();

  std::sort(events.begin(), events.end(), comp);
  vector<segment2D> reported;
  for (int i = 0; i < events.size(); i++) {
    sweepEvent curr = events.at(i);

    if (curr.type == 0) {
      as.insert(std::pair<int, segment2D>(curr.s.start.y, curr.s));
    }
    else if (curr.type == 1) {
        std::pair<multimap<int, segment2D>::iterator, multimap<int, segment2D>::iterator> ret;
        ret = as.equal_range(curr.s.end.y);
        
        for (multimap<int, segment2D>::iterator it = ret.first; it != ret.second; ++it) {
            if (curr.s.end.x == it->second.end.x) {
                as.erase(it);
                break;
            }
        }
        
        
    }
    else {
      int y1, y2;

      if (curr.s.start.y > curr.s.end.y) {
        y1 = curr.s.end.y;
        y2 = curr.s.start.y;
      } else {
        y1 = curr.s.start.y;
        y2 = curr.s.end.y;
      }
      
      reported = rangeSearch(y1, y2);
      addPoints(reported, curr.s.start.x);
      reported.clear();
    }
  }
}

void addPoints(vector<segment2D> s, int x) {
  int y;
  for (int i = 0; i < s.size(); i++) {
    y = s[i].start.y;
    point2D p;
    p.x = x;
    p.y = y;
    intpoints.push_back(p);
  }
}

vector<segment2D> rangeSearch(int y1, int y2) {
  vector<segment2D> inRange;
  multimap<int, segment2D>::iterator it1, it2;
  it1 = as.lower_bound(y1);
  it2 = as.upper_bound(y2);
  for (multimap<int, segment2D>::iterator it = it1; it != it2; ++it) {
    inRange.push_back(it->second);
  }

  return inRange;
}

bool comp(const sweepEvent &a, const sweepEvent &b) { return a.x < b.x; }

void getEvents() {
  for (int i = 0; i < n; i++) {
    if (segments[i].start.y == segments[i].end.y) {
      sweepEvent event1;
      event1.x = segments[i].start.x;
      event1.type = 0;
      event1.s = segments[i];

      sweepEvent event2;
      event2.x = segments[i].end.x;
      event2.type = 1;
      event2.s = segments[i];

      events.push_back(event1);
      events.push_back(event2);
    } else {
      sweepEvent event;
      event.x = segments[i].start.x;
      event.type = 2;
      event.s = segments[i];

      events.push_back(event);
    }
  }
}


/* ************************************************** */
//fills global variable "segments" with n segments 
void initialize_segments() {

  switch (init_case)  {
      
    case 0: 
      initialize_segments_random(); 
      break;
      
    case 1: 
      initialize_segments_horizontal(); 
      break;
          
    case 2:
      initialize_segments_overlap();
      break;
          
    case 3:
      initialize_segments_grid();
      break;
          
    case 4:
      initialize_segments_vertical();
      break;
      
    default: 
      initialize_segments_random(); 
    }

  init_case = (init_case+1) % NB_TEST_CASES;
  return; 
}


//Marcus's test cases
/* ************************************************** */
/* ************************************************** */

void initialize_segments_overlap_vertical(){
  
  intpoints.clear();
  segments.clear();
  point2D a,b;
  segment2D s;

  for (int i = 1; i < 4; i++) {
    a.x = 0;
    b.x = WINDOWSIZE;
    a.y = 1 + i * WINDOWSIZE / 4;
    b.y = a.y;
    
    s.start = a;
    s.end = b;
    segments.push_back(s);
  }

  double k = 0.25;
    
  for (int i = 1; i < 4; i++) {
    for (int j = 1; j < 3; j++) {
      a.x = WINDOWSIZE*k;
      b.x = WINDOWSIZE*k;
      a.y = 0;
      b.y = WINDOWSIZE;
      
      s.start = a;
      s.end = b;
      segments.push_back(s);
    }

  k += 0.25;

  }
    
}

void initialize_segments_overlap_horizontal(){
  
  intpoints.clear();
  segments.clear();
  point2D a,b;
  segment2D s;

  int k = 0;

  for (int i = 0; i < 10; i++) {
    a.x = 0;
    b.x = k + 50;
    a.y = WINDOWSIZE/2;
    b.y = a.y;
    
    s.start = a;
    s.end = b;
    segments.push_back(s);

    k += 50;
  }
    
  for (int i = 1; i < 4; i++) {
    a.x = 1 + i * WINDOWSIZE / 4;
    b.x = a.x;
    a.y = 1;
    b.y = WINDOWSIZE - 1;
    
    s.start = a;
    s.end = b;
    segments.push_back(s);
  }
   
}

void initialize_segments_lattice(){
  
  intpoints.clear();
  segments.clear();
  point2D a,b;
  segment2D s;

  int k = 10;

  for (int i = 0; i < 49; i++) { //Horizontal
    a.x = 0;
    b.x = WINDOWSIZE;
    a.y = k;
    b.y = a.y;
        
    s.start = a;
    s.end = b;
    segments.push_back(s);

    k += 10;

  }

  k = 10;

  for (int i = 0; i < 49; i++) { //Vertical
    a.x = k;
    b.x = k;
    a.y = 0;
    b.y = WINDOWSIZE;
        
    s.start = a;
    s.end = b;
    segments.push_back(s);

    k += 10;

  }

    
}

void initialize_segments_rolling(){
  
  intpoints.clear();
  segments.clear();
  point2D a,b;
  segment2D s;

  int k = 0;

  for (int i = 0; i < 9; i++) { //Horizontal
    a.x = k;
    b.x = a.x + 100;
    a.y = WINDOWSIZE / 2;
    b.y = a.y;

    s.start = a;
    s.end = b;
    segments.push_back(s);

    k += 50;

  }

  k = 50;

  for (int i = 0; i < 9; i++) { //Vertical
    a.x = k;
    b.x = k;
    a.y = 0;
    b.y = WINDOWSIZE;
        
    s.start = a;
    s.end = b;
    segments.push_back(s);

    k += 50;

  }
   
}

/* ************************************************** */
//fills global variable "segments" with n horizontal segments 
void initialize_segments_horizontal() {
    
  int i; 
  point2D a,b;
  segment2D s; 

  //clear the vector
  segments.clear();
  intpoints.clear();

  //a long horizontal segment 
  a.x = 1; 
  a.y = WINDOWSIZE/2; 
  b.x = WINDOWSIZE - 10; 
  b.y = a.y; 

  s.start = a; s.end = b; 
  segments.push_back(s);  

  //n-1 vertical segments 
  for (i=0; i<n-1; i++) {
    
    a.x = i*WINDOWSIZE/n; 
    a.y = WINDOWSIZE/2 - random() % ((int)(.4*WINDOWSIZE)); 
    b.x = a.x; 
    b.y = WINDOWSIZE/2 + random() % ((int)(.4*WINDOWSIZE)); 
    s.start = a; s.end = b; 
    segments.push_back(s); 
  }

}


/* ****************************** */
//fills global variable "segments" with n random segments 
void initialize_segments_random() {
  
  //clear the vector 
  segments.clear(); 
  intpoints.clear();

  int i; 
  point2D a, b; 
  segment2D s; 
  for (i=0; i<n; i++) {
    if (random()%2 == 0) {
      //horizontal segment
      a.x = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE)); 
      a.y =  (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE));
      b.y = a.y; 
      b.x = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE)); 

      if (a.x < b.x) {
	s.start = a; s.end = b; 
      } else {
	s.start = b; s.end = a; 
      } 
 
    } else {
      //vertical segment 
      a.x = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE)); 
      b.x = a.x; 
      a.y = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE)); 
      b.y = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE)); 

      if (a.y < b.y) {
	s.start = a; s.end = b; 
      } else {
	s.start = b; s.end = a; 
      }
    }

    //insert the segment in the array of segments 
    segments.push_back (s); 
  } //for i
}

//fills global variable "segments" with 1 vertical segment and n-1
//horizontal segments
void initialize_segments_vertical() {
    
    int i;
    point2D a,b;
    segment2D s;
    //clear the vector
    segments.clear();
    intpoints.clear();
    
    //a long vertical segment
    a.y = 1;
    a.x = WINDOWSIZE/2;
    b.y = WINDOWSIZE - 10;
    b.x = a.x;
    
    s.start = a; s.end = b;
    segments.push_back(s);
    
    //n-1 horizontal segments
    for (i=0; i<n-1; i++) {
        
        a.y = i*WINDOWSIZE/n;
        a.x = WINDOWSIZE/2 - random() % ((int)(.4*WINDOWSIZE));
        b.y = a.y;
        b.x = WINDOWSIZE/2 + random() % ((int)(.4*WINDOWSIZE));
        
        if (a.x < b.x) {
            s.start = a; s.end = b;
        } else {
            s.start = b; s.end = a;
        }
        
        segments.push_back(s);
    }
    
}

//fills the global variable "segments" with an n x n grid
void initialize_segments_grid() {
    int i;
    point2D a,b;
    segment2D s;
    
    //clear the vector
    segments.clear();
    
    //n-1 vertical segments
    for (i=0; i < n; i++) {
        if (i % 2 == 0) {
            a.x = (2 * i) % (WINDOWSIZE / 2) + (WINDOWSIZE / 4);
            a.y = (WINDOWSIZE / 4);
            b.x = a.x;
            b.y = (WINDOWSIZE / 4) + (WINDOWSIZE / 2);
        }
        else {
            a.y = (2 * i) % (WINDOWSIZE / 2) + WINDOWSIZE / 4;
            a.x = (WINDOWSIZE / 4);
            b.y = a.y;
            b.x = (WINDOWSIZE / 4) + (WINDOWSIZE / 2);
        }
        s.start = a; s.end = b;
        segments.push_back(s);
    }
}

//fills the global variable "segments" with a 3x3 grid of segments that meet
//at their endpoints
void initialize_segments_overlap() {
    int i;
    point2D a,b;
    segment2D s;

    intpoints.clear();
    segments.clear();
    
    //n-1 vertical segments
    for (i=0; i < n; i++) {
        if (i % 6 == 0) {
            a.x = WINDOWSIZE/4;
            a.y = WINDOWSIZE/4;
            b.x = a.x;
            b.y = 3*WINDOWSIZE/4;
        }
        else if (i % 6 == 1) {
            a.x = WINDOWSIZE/2;
            a.y = WINDOWSIZE/4;
            b.x = a.x;
            b.y = 3*WINDOWSIZE/4;
        }
        else if (i % 6 == 2) {
            a.x = 3*WINDOWSIZE/4;
            a.y = WINDOWSIZE/4;
            b.x = a.x;
            b.y = 3*WINDOWSIZE/4;
        }
        else if (i % 6 == 3) {
            a.y = WINDOWSIZE/4;
            a.x = WINDOWSIZE/4;
            b.y = a.y;
            b.x = 3*WINDOWSIZE/4;
        }
        else if (i % 6 == 4) {
            a.y = WINDOWSIZE/2;
            a.x = WINDOWSIZE/4;
            b.y = a.y;
            b.x = 3*WINDOWSIZE/4;
        }
        else {
            a.y = 3*WINDOWSIZE/4;
            a.x = WINDOWSIZE/4;
            b.y = a.y;
            b.x = 3*WINDOWSIZE/4;
        }
        
        s.start = a; s.end = b; 
        segments.push_back(s); 
    }
}



/* ************************************************** */
void print_segments() {

  for (int i=0; i<segments.size(); i++) {
    printf("segment %d: [(%d,%d), (%d,%d)]\n",
	   i, segments[i].start.x, segments[i].start.y, segments[i].end.x, segments[i].end.y);
  }
}




/* ****************************** */
int main(int argc, char** argv) {

  //read number of points from user
  if (argc!=2) {
    printf("usage: viewPoints <nbPoints>\n");
    exit(1); 
  }
  n = atoi(argv[1]); 
  printf("you entered n=%d\n", n);
  assert(n >0); 

  initialize_segments_random();

  Rtimer rt1;
  rt_start(rt1);

  //compute something here 
  intersection();

  printf("Number of Intersections: %lu\n", intpoints.size());


  rt_stop(rt1);

  //print the timing 
  char buf [1024];
  rt_sprint(buf,rt1);
  printf("run time:  %s\n\n", buf);
  fflush(stdout);

 
  

  /* initialize GLUT  */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(WINDOWSIZE, WINDOWSIZE);
  glutInitWindowPosition(100,100);
  glutCreateWindow(argv[0]);

  /* register callback functions */
  glutDisplayFunc(display); 
  glutKeyboardFunc(keypress);
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
/* draw the segments stored in global variable segments */
void draw_segments(){

  //set color 
  glColor3fv(white);
  glLineWidth(1);   
  
  int i;
  for (i=0; i<segments.size(); i++) {
    glBegin(GL_LINES);
    glVertex2f(segments[i].start.x, segments[i].start.y); 
    glVertex2f(segments[i].end.x, segments[i].end.y);
    glEnd();
  }
}

//draw the sweep line 
void draw_sweep_line() {

  //sweep line color 
  glColor3fv(green); 

  //the current position of sweep line is sweep_line_x; assume it's a
  //segment from y=0 to y=windowsize;
  glBegin(GL_LINES); 
  glVertex2f(sweep_line_x, 0); 
  glVertex2f(sweep_line_x, WINDOWSIZE); 
  glEnd();
}

//draw a segment with current color 
void draw_segment(segment2D s, GLfloat width) {
  
  glLineWidth(width);
  glBegin(GL_LINES);
  glVertex2f(s.start.x, s.start.y); 
  glVertex2f(s.end.x, s.end.y);
  glEnd();
}

//draw all the elements in the active structure 
void draw_active_structure() {

  glColor3fv(red);

  for (int i = 0; i < segments.size(); i++) {

    segment2D curr = segments[i];
    int x1, x2;

    if (curr.start.x >= curr.end.x) {
      x1 = curr.end.x;
      x2 = curr.start.x;
    } else {
      x1 = curr.start.x;
      x2 = curr.end.x;
    }

    if (x1 <= sweep_line_x && sweep_line_x <= x2) {
      draw_segment(curr, 3);
    }
  }
}


//draw all the elements in intpoints 
void draw_intersection_points() {

  glColor3fv(cyan);
  const int R = 2;

  for (int i = 0; i < intpoints.size(); i++) {
    
    if (sweep_line_x >= intpoints[i].x) {
      
      glBegin(GL_POLYGON);
      glVertex2f(intpoints[i].x -R, intpoints[i].y-R);
      glVertex2f(intpoints[i].x +R, intpoints[i].y-R);
      glVertex2f(intpoints[i].x +R, intpoints[i].y+R);
      glVertex2f(intpoints[i].x -R, intpoints[i].y+R);
      glEnd();
    }
  }
}



/* ****************************** */
void display(void) {

  glClear(GL_COLOR_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); //clear the matrix


  /* The default GL window is [-1,1]x[-1,1] with the origin in the
     center. 
     
     The points are in the range (0,0) to (WINSIZE,WINSIZE), so they
     need to be mapped to [-1,1]x [-1,1]x */
  glScalef(2.0/WINDOWSIZE, 2.0/WINDOWSIZE, 1.0);  
  glTranslatef(-WINDOWSIZE/2, -WINDOWSIZE/2, 0); 


  draw_active_structure(); 
  draw_segments();
  draw_intersection_points(); 
  draw_sweep_line(); 


  /* execute the drawing commands */
  glFlush();
}



/* ****************************** */
void keypress(unsigned char key, int x, int y) {
  switch(key) {
  case 'q':
    exit(0);
    break;

  case 'i': 
    initialize_segments(); 
    glutPostRedisplay();
    break;
  

          
  case 'r':
    Rtimer rt1;
    rt_start(rt1);

    intersection();
          

    rt_stop(rt1);
          
    //print the timing
    printf("Number of Intersections: %lu\n", intpoints.size());
    char buf [1024];
    rt_sprint(buf,rt1);
    printf("run time:  %s\n\n", buf);
    fflush(stdout);
    
    sweep_line_x = 0;
    glutDisplayFunc(display);
    
    break;
          
      case 'c':
          intpoints.clear();
          glClear(GL_COLOR_BUFFER_BIT);
          glMatrixMode(GL_MODELVIEW);
          glLoadIdentity();
          glutPostRedisplay();
          
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
  
  /* LT: I used this to slow things down in a controlled way; probably not
     necessary for this assignment..*/
  //static int lastFrameTime=0;  
  //note: a static variable, remembered from one call to the next
  //int now, elapsed_ms; 
  
  //now = glutGet (GLUT_ELAPSED_TIME); 
  //elapsed_ms = now - lastFrameTime; 
  //lastFrameTime=now; 
  
  sweep_line_x++; 

  
  glutPostRedisplay(); 

}
