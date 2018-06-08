CC = gcc
CFLAGS = -g -O3 -lm -lgsl -lgslcblas -w

simSFS: 
	$(CC) -o simSFS *.c $(CFLAGS)

clean: 
		rm simSFS *.o