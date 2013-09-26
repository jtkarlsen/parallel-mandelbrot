CC = mpicc -O3

all: libgraphicsscreen.a libStopWatch.a

libgraphicsscreen.a: graphicsScreen.o
	ar cru libgraphicsscreen.a graphicsScreen.o

graphicsScreen.o: graphicsScreen.c
	$(CC) $(INC) -c graphicsScreen.c

libStopWatch.a: StopWatch.o
	ar cru libStopWatch.a StopWatch.o

StopWatch.o: StopWatch.c
	$(CC) $(INC) -c StopWatch.c

clean:
	rm -f *.o *~ *.a
