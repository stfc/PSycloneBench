
# The targets that this Makefile supports
EXECS = nemolite2d

all: $(EXECS)

COMMON_MODULES = continuity_kern.o

# Normal targets
.PHONY: nemolite2d
nemolite2d: 
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES}" nemolite2d.exe

${API_LIB}: ${API_DIR}/*.?90
	${MAKE} -C ${API_DIR} F90="${F90}" F90FLAGS="${F90FLAGS} ${OMPFLAGS}" AR="${AR}" ARFLAGS="${ARFLAGS}" API_LIB="gocean_api.a"

nemolite2d.o: $(COMMON_MODULES)

# Interdependencies between modules, alphabetical order

# Generic rules

%.exe: $(MODULE_LIST)
	$(CC) -o $@ $(MODULE_LIST) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) ${OMPFLAGS} -c $<

clean: 
	rm -f *.o *~

allclean: clean
	rm -f *.exe fparser.log
	rm -rf *_wpl_dir/
