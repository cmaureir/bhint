
# +--------------------------------------------+
# + Makefile for integrators                   +
# + U. Loeckmann                               +
# +--------------------------------------------+

LIBS=

CFLAGS=-Wall -O3 -pg -fopenmp
CC=g++
NVCC=nvcc
GET=co
OBJS=bhi_int.o  bhi_io.o  bhi_kepler.o  bhi_timestep.o  bhi_util.o  bhi_vector.o
DEPEND= makedepend $(CFLAGS)
LDFLAGS=-L . -L /opt/SUNWhpc/lib $(G2C) $(LIBS) -lm -pg -fopenmp -Wall
CUFLAGS= -L /usr/local/cuda/lib64 -lcudart -lcuda


all: bhint

bhi_main.o: bhi_main.c
	$(CC) $(CFLAGS) -c bhi_main.c

clean:
	rm -f *.o bhint

cleanall: clean
	make -C clean

bhi_int.o: bhi_int.cu
	$(NVCC) -Xptxas -v -c $^ -arch sm_13

%.o: %.c bhi.h bhi_config.h
	$(CC) $(CFLAGS) -c $*.c

# To make an executable

bhint: $(OBJS) bhi_main.o
	$(CC) -o bhint $(OBJS) bhi_main.o $(LDFLAGS) $(CUFLAGS)
