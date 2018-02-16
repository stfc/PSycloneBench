CLFORTRAN_PATH = ${HOME}/Projects/clfortran
CLFORTRAN_INCLUDE = -I${CLFORTRAN_PATH}

COMMON_MODULES = opencl_utils_mod.o nemolite2d.o
nemolite2d.exe: ${COMMON_MODULES}
	${F90} -o $@ ${COMMON_MODULES} ${OPENCL_LIBS} $(LDFLAGS)

%.o: %.f90
	${F90} ${F90FLAGS} -c ${CLFORTRAN_INCLUDE} $<

clean:
	rm -f *.o *~

allclean: clean
	rm -f *.exe *.mod
