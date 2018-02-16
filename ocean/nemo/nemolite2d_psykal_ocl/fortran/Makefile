CLFORTRAN_PATH = ${HOME}/Projects/clfortran
CLFORTRAN_INCLUDE = -I${CLFORTRAN_PATH}

nemolite2d.exe: nemolite2d.o
	${F90} -o $@ nemolite2d.o ${CLFORTRAN_PATH}/clfortran.o ${OPENCL_LIBS} $(LDFLAGS)

%.o: %.f90
	${F90} ${F90FLAGS} -c ${CLFORTRAN_INCLUDE} $<

clean:
	rm -f *.o *~

allclean: clean
	rm *.exe
