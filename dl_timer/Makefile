# The name of the library/archive file we will build
DLT_LIB = dl_timer_lib.a
# Location of the source files
SRC_DIR = src
# Location of the test code
TEST_DIR = test

all:
	@echo "Possible make targets are:"
	@echo "   Build: sm_lib (OpenMP support), dm_lib (MPI support)"
	@echo "   Test: sm_test, dm_test"

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

# The directory 'test' does actually exist but this target does not
# create or update it - therefore mark it as phony. By default we do the
# non-MPI tests.
.PHONY: test
test: sm_test

.PHONY: sm_test
sm_test: sm_lib
	${MAKE} --directory=test sm_test

.PHONY: dm_test
dm_test: dm_lib
	${MAKE} --directory=test dm_test

.PHONY: testclean
testclean:
	${MAKE} --directory=${TEST_DIR} clean

.PHONY: clean
clean: testclean
	${MAKE} --directory=${SRC_DIR} clean

.PHONY: allclean
allclean:
	${MAKE} --directory=${SRC_DIR} allclean
	${MAKE} --directory=${TEST_DIR} allclean
	rm -f ${DLT_LIB} *~
