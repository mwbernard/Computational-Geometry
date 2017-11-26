PLATFORM = $(shell uname)


## Compilation flags
##comment out one or the other 
##debugging 
CFLAGS = -g 
##release
#CFLAGS = -O3 -DNDEBUG
LDFLAGS=

CFLAGS+= -Wall

ifeq ($(PLATFORM),Darwin)
## Mac OS X
CFLAGS += -m64 -isystem/usr/local/include  -Wno-deprecated 
LDFLAGS+= -m64 -lc -framework AGL -framework OpenGL -framework GLUT -framework Foundation
else
## Linux
CFLAGS += -m64
INCLUDEPATH  = -I/usr/include/GL/ 
LIBPATH = -L/usr/lib64 -L/usr/X11R6/lib
LDFLAGS+=  -lGL -lglut -lrt -lGLU -lX11 -lm  -lXmu -lXext -lXi
endif


CC = gcc $(CFLAGS) $(INCLUDEPATH)


PROGS = viewPoints

default: $(PROGS)

viewPoints: viewPoints.o geom.o rtimer.o  kdtree.o
	$(CC) -o $@ viewPoints.o geom.o rtimer.o  kdtree.o $(LDFLAGS)

viewPoints.o: viewPoints.c kdtree.h rtimer.h
	$(CC) -c $(INCLUDEPATH)  viewPoints.c  -o $@

geom.o: geom.c geom.h 
	$(CC) -c $(INCLUDEPATH)  geom.c -o $@

rtimer.o: rtimer.h rtimer.c
	$(CC) -c $(INCLUDEPATH) rtimer.c -o $@

kdtree.o: kdtree.c kdtree.h geom.h
	$(CC) -c $(INCLUDEPATH)  kdtree.c -o $@

clean::	
	rm *.o
	rm viewPoints


