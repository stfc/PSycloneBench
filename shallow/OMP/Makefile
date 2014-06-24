# Makefile.include is a symbolic link to the desired
# Makefile.include.<compiler> file.
include Makefile.include

all: shallow_omp.exe

clean:
	rm *.exe
	rm *~
	rm *.o

OBJECTS = wtime.o shallow_base_openmp_v3.o

shallow_omp.exe: ${OBJECTS}
	${CC} ${CFLAGS} ${OMPFLAGS} -o $@ ${OBJECTS} ${LDFLAGS}

%.o: %.c
	${CC} ${CFLAGS} ${OMPFLAGS} -c $<
