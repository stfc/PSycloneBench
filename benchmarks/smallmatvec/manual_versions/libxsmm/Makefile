FC=ifort
FFLAGS=-qopenmp -Ofast
LDFLAGS=

# MKLFLAGS="-L${MKLROOT}/lib/intel64 -mkl"
XSMM_INCLUDE=-I${HOME}/libxsmm/libxsmm-1.15/include
XSMM_LINK=-L/home/rford/libxsmm/libxsmm-1.15/lib -lxsmmf -lxsmm -lxsmmnoblas

all : kdriver_orig kdriver_xsmm

# Dependencies
kdriver_orig.o: matrix_vector_kernel_mod.o constants_mod.o
matrix_vector_kernel_mod.o: constants_mod.o

kdriver_xsmm.o: matrix_vector_kernel_xsmm_mod.o constants_mod.o
matrix_vector_kernel_xsmm_mod.o: constants_mod.o

# Compile object targets
%.o: %.F90
	$(FC) $(FFLAGS) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

kdriver_orig: kdriver_orig.o matrix_vector_kernel_mod.o constants_mod.o
	$(FC) $(FFLAGS) -xAVX2 -o $@ $^ $(LDFLAGS)

kdriver_orig.o: kdriver.f90
	$(FC) $(FFLAGS) -c $< -o kdriver_orig.o

kdriver_xsmm: kdriver_xsmm.o matrix_vector_kernel_xsmm_mod.o constants_mod.o
	$(FC) $(FFLAGS) -xAVX2 -o $@ $^ $(LDFLAGS) ${XSMM_LINK}

kdriver_xsmm.o: kdriver.f90
	$(FC) $(FFLAGS) -c $< -o kdriver_xsmm.o

matrix_vector_kernel_xsmm_mod.o: matrix_vector_kernel_xsmm_mod.f90
	$(FC) $(FFLAGS) -c $< -o matrix_vector_kernel_xsmm_mod.o ${XSMM_INCLUDE}

run:
	@for v in $(VERSIONS); do \
		./kdriver; \
	done

versions:
	@echo $(VERSIONS)

clean:
	rm -f *.o *.mod *.optrpt *.s kdriver_orig kdriver_xsmm

allclean: clean

