# Makefile for the C, OpenMP version of Shallow. Compiler and associated
# flags are picked up from environment variables. See e.g.
# ../../../compiler_setup/gnu.sh

all: shallow_omp.exe

clean:
	rm -f *~
	rm -f *.o

allclean: clean
	rm -f *.exe

OBJECTS = wtime.o shallow_base_openmp_v3.o

LDFLAGS += -lm

shallow_omp.exe: ${OBJECTS}
	${CC} ${CFLAGS} ${OMPFLAGS} -o $@ ${OBJECTS} ${LDFLAGS}

%.o: %.c
	${CC} ${CFLAGS} ${OMPFLAGS} -c $<
