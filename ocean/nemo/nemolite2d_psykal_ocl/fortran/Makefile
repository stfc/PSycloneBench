CLFORTRAN_PATH = ${HOME}/Projects/clfortran
CLFORTRAN_INCLUDE = -I${CLFORTRAN_PATH}

# Location of the code defining the GOcean API
API_VERSION=1.0
API_DIR = ../../../api_v${API_VERSION}
API_LIB = ${API_DIR}/gocean_api.a

# Location of the dl_timer library
TIMER_LIB = ${HOME}/Projects/dl_timer/dl_timer_lib.a
TIMER_INC = ${HOME}/Projects/dl_timer/src

COMMON_MODULES = opencl_utils_mod.o kernel_args_mod.o gocean2d_io_mod.o \
                 model_mod.o initialisation_mod.o nemolite2d.o ${API_LIB}

nemolite2d.exe: ${API_LIB} ${COMMON_MODULES}
	${F90} -o $@ ${COMMON_MODULES} ${OPENCL_LIBS} ${TIMER_LIB} $(LDFLAGS)

%.o: %.f90
	${F90} ${F90FLAGS} -c ${CLFORTRAN_INCLUDE} -I${API_DIR} -I${TIMER_INC} $<

${API_LIB}: ${API_DIR}/*.?90
	${MAKE} -C ${API_DIR} F90="${F90}" F90FLAGS="${F90FLAGS}" \
            OMPFLAGS="${OMPFLAGS}" AR="${AR}" ARFLAGS="${ARFLAGS}" \
            API_LIB="gocean_api.a"

clean:
	rm -f *.o *~

allclean: clean
	rm -f *.exe *.mod
