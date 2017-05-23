CC = gcc
CFLAGS = -g -O2 -lm -lgsl -lgslcblas -w

simSFS: simSFS.v1.0.o
	$(CC) -o simSFS simSFS.v1.0.o $(CFLAGS)

clean: 
		rm simSFS *.o
