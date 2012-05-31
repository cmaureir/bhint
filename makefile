
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

CFLAGS=-Wall -O3 -I /opt/SUNWhpc/include
CC=gcc
CC_FORTRAN=gcc
GET=co
OBJS=bhi_int.o  bhi_io.o  bhi_kepler.o  bhi_sse.o  bhi_timestep.o  bhi_util.o  bhi_vector.o

DEPEND= makedepend $(CFLAGS)

LDFLAGS=-L . -L /opt/SUNWhpc/lib $(G2C) $(LIBS) -lm

default: bhint

bhi_main.o: bhi_main.c FORCE
	$(CC) $(CFLAGS) -c bhi_main.c


FORCE:
	
clean:
	rm -f *.o bhint

cleanall: clean
	make -C sse clean

%.o: %.c bhi.h bhi_config.h
	$(CC) $(CFLAGS) -c $*.c

libsse.a:
	@case "$(LIBS)" in *-lsse* ) make -C sse ;; esac

# To make an executable

bhint: $(OBJS) bhi_main.o libsse.a
	$(CC_FORTRAN) -o bhint $(OBJS) bhi_main.o $(LDFLAGS2) $(LDFLAGS)

