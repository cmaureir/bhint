
# +--------------------------------------------+
# + Makefile for integrators                   +
# + U. Loeckmann                               +
# +--------------------------------------------+

# additional libraries (stellar evolution, GRAPE etc.)
LIBS=
#LIBS=-lsse -lgfortran
#LIBS=-lg2c -lsse
#LIBS=-lg2c -lsse -lg6lx
#LIBS=-lg2c -lsse -lg6a_new

CFLAGS=-Wall -O3 -pg -fopenmp
CC=gcc
CC_FORTRAN=gcc
GET=co
OBJS=bhi_int.o  bhi_io.o  bhi_kepler.o  bhi_sse.o  bhi_timestep.o  bhi_util.o  bhi_vector.o
DEPEND= makedepend $(CFLAGS)
LDFLAGS=-L . -L /opt/SUNWhpc/lib $(G2C) $(LIBS) -lm -pg -fopenmp


all: bhint

bhi_main.o: bhi_main.c
	$(CC) $(CFLAGS) -c bhi_main.c

clean:
	rm -f *.o bhint

cleanall: clean
	make -C sse clean

%.o: %.c bhi.h bhi_config.h
	$(CC) $(CFLAGS) -c $*.c

libsse.a:
	@case "$(LIBS)" in *-lsse* ) make -C ../original/sse ;; esac

# To make an executable

bhint: $(OBJS) bhi_main.o libsse.a
	$(CC_FORTRAN) -o bhint $(OBJS) bhi_main.o $(LDFLAGS2) $(LDFLAGS)

