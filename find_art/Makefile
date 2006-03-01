PROGS=find_art
OBJS=$(PROGS:=.o) \
	../minc_vector_io.o \
	../perf_util.o \
	callbacks.o \
	interface.o \
	gtk_gl.o


CC = gcc

OPTIONS = -g -O3 -Wall -Wunused -Wmissing-prototypes -Wmissing-declarations

INCLUDES = -I/usr/local/mni/include \
`gsl-config --cflags` \
`gtk-config --cflags`

LDINCLUDES = -L/usr/local/mni/lib \
`gsl-config --libs` \
`gtk-config --libs`

LDLIBS = -lgtkgl -lGL -lGLU -lpthread -lvolume_io -lminc -lnetcdf -lm

LDOPTS = $(LDINCLUDES) $(LDLIBS)
CFLAGS = $(OPTIONS) $(INCLUDES)


all: $(PROGS) 

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGS): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) -o $@ $(LDOPTS)
        
clean:
	rm -f *.o *~ $(PROGS)
