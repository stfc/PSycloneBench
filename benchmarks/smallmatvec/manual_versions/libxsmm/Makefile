FC=ifort
FFLAGS=-qopenmp -Ofast
LDFLAGS=

# MKLFLAGS="-L${MKLROOT}/lib/intel64 -mkl"
XSMM_INCLUDE=-I/lustre/scafellpike/local/HCP055/pxs02/rxf70-pxs02/libxsmm-1.15/include
XSMM_LINK=-L/lustre/scafellpike/local/HCP055/pxs02/rxf70-pxs02/libxsmm-1.15/lib -lxsmmf -lxsmm -lxsmmnoblas

all : kdriver_orig kdriver_xsmm kdriver_reorder

# Dependencies
kdriver_orig.o: matrix_vector_kernel_mod.o constants_mod.o
matrix_vector_kernel_mod.o: constants_mod.o

kdriver_xsmm.o: matrix_vector_kernel_xsmm_mod.o constants_mod.o
matrix_vector_kernel_xsmm_mod.o: constants_mod.o

kdriver_reorder.o: matrix_vector_kernel_reorder_mod.o constants_mod.o
matrix_vector_kernel_reorder_mod.o: constants_mod.o

# Compile object targets
%.o: %.F90
	$(FC) $(FFLAGS) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

kdriver_orig: kdriver_orig.o matrix_vector_kernel_mod.o constants_mod.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

kdriver_orig.o: kdriver.f90
	$(FC) $(FFLAGS) -c $< -o kdriver_orig.o

kdriver_xsmm: kdriver_xsmm.o matrix_vector_kernel_xsmm_mod.o constants_mod.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) ${XSMM_LINK}

kdriver_xsmm.o: kdriver.f90
	$(FC) $(FFLAGS) -c $< -o kdriver_xsmm.o

matrix_vector_kernel_xsmm_mod.o: matrix_vector_kernel_xsmm_mod.f90
	$(FC) $(FFLAGS) -c $< -o matrix_vector_kernel_xsmm_mod.o ${XSMM_INCLUDE}

kdriver_reorder: kdriver_reorder.o matrix_vector_kernel_reorder_mod.o constants_mod.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

kdriver_reorder.o: kdriver.f90
	$(FC) $(FFLAGS) -c $< -o kdriver_reorder.o

matrix_vector_kernel_reorder_mod.o: matrix_vector_kernel_reorder_mod.f90
	$(FC) $(FFLAGS) -c $< -o matrix_vector_kernel_reorder_mod.o ${XSMM_INCLUDE}

run:
	@for v in $(VERSIONS); do \
		./kdriver; \
	done

versions:
	@echo $(VERSIONS)

clean:
	rm -f *.o *.mod *.optrpt *.s kdriver_orig kdriver_xsmm kdriver_reorder

allclean: clean

