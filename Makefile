PROGS = volperf
HEADERS = perf_util.h gamma_fit.h minc_vector_io.h 
OBJS = $(PROGS:=.o) $(HEADERS:.h=.o)

OPTIONS = -O3
INCLUDES = -I/usr/local/include -I/usr/local/mni/include
CFLAGS = $(OPTIONS) $(INCLUDES) `gsl-config --cflags`

LDINCLUDES = -L/usr/local/lib32 -L/usr/local/mni/lib32
LDLIBS = -lvolume_io -lminc -lnetcdf -lm
LDOPTS = $(LDINCLUDES) $(LDLIBS) `gsl-config --libs`

all: $(PROGS) 

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGS): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) -o $@ $(LDOPTS)

clean:
	rm -f *.o *~ $(PROGS)
