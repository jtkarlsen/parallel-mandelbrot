CC = mpicc
# Remove -DGRAPHICS to compile with no graphic support (e.g. for measurements)
CFLAGS = -O2 -Wall -g -DGRAPHICS
RRDIR = ../lib

# Adapt this to your distribution if necessary
INCSDL =
LIBSDL = -L/usr/X11R6/lib64 


INC = $(INCSDL) -I$(RRDIR) 
LIB = $(LIBSDL) -L$(RRDIR) -lX11 -lgraphicsscreen -lSDL -lStopWatch -lpthread -lm

all: RoadMap

RoadMap: RoadMap.c
	$(CC) $(CFLAGS) RoadMap.c -o RoadMap $(INC) $(LIB)

clean:
	rm -f RoadMap RoadMapOMP RoadMapMPI *~
