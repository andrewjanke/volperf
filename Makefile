PROGS = volperf
HEADERS = perf_util.o minc_vector_io.o 
OBJS = $(PROGS:=.o) $(HEADERS:.h=.o)

CC=cc

OPTIONS = -g3 -fullwarn -O3
INCLUDES = -I/usr/local/include
CFLAGS = $(OPTIONS) $(INCLUDES) `gsl-config --cflags`

LDINCLUDES = -L/usr/local/lib32
LDLIBS = -lvolume_io -lminc -lnetcdf -lm
LDOPTS = $(LDINCLUDES) $(LDLIBS) `gsl-config --libs`

all: $(PROGS) 

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGS): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) -o $@ $(LDOPTS)

clean:
	rm -f *.o *~ $(PROGS)
