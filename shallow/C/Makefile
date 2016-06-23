# Makefile.include is a symbolic link to the desired
# Makefile.include.<compiler> file.
include Makefile.include

all: shallow_base

clean:
	rm -f *.exe
	rm -f *~
	rm -f *.o

OBJECTS = shallow_base.o

shallow_base: ${OBJECTS}
	${CC} ${CFLAGS} ${OMPFLAGS} -o $@ ${OBJECTS} ${LDFLAGS}

%.o: %.c
	${CC} ${CFLAGS} ${OMPFLAGS} -c $<
