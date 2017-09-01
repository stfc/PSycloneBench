# The name of the library/archive file we will build
DLT_LIB = dl_timer_lib.a
# Location of the source files
SRC_DIR = src

all:
	@echo "Possible make targets are:"
	@echo "   Build: sm_lib (OpenMP support), dm_lib (MPI support)"

# Library with shared memory (OpenMP) support
.PHONY: sm_lib
sm_lib: ${SRC_DIR}/*.?90
	${MAKE} --directory=${SRC_DIR} LIB_NAME=${DLT_LIB} sm_build
	mv ${SRC_DIR}/${DLT_LIB} .

# Library with distributed memory (MPI) support
.PHONY: dm_lib
dm_lib: ${SRC_DIR}/*.?90
	${MAKE} --directory=${SRC_DIR} LIB_NAME=${DLT_LIB} dm_build
	mv ${SRC_DIR}/${DLT_LIB} .

.PHONY: clean
clean:
	${MAKE} --directory=${SRC_DIR} clean

.PHONY: allclean
allclean:
	${MAKE} --directory=${SRC_DIR} allclean
	rm -f ${DLT_LIB} *~
