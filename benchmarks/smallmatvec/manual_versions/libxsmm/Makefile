FC ?= ifort
FFLAGS ?= -qopenmp -Ofast
LDFLAGS ?=

XSMM_DIR ?= ${HOME}/libxsmm-1.15
XSMM_INCLUDE ?= -I${XSMM_DIR}/include
XSMM_LINK ?= -L${XSMM_DIR}/lib -lxsmmf -lxsmm -lxsmmnoblas

all : kdriver_orig kdriver_xsmm kdriver_reorder

# Dependencies
kdriver_orig.o: matrix_vector_kernel_mod.o constants_mod.o
matrix_vector_kernel_mod.o: constants_mod.o

kdriver_xsmm.o: matrix_vector_kernel_xsmm_mod.o constants_mod.o
matrix_vector_kernel_xsmm_mod.o: constants_mod.o

kdriver_reorder.o: matrix_vector_kernel_reorder_mod.o constants_mod.o
matrix_vector_kernel_reorder_mod.o: constants_mod.o

# Compile object targets
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
	$(FC) $(FFLAGS) -c $< -o matrix_vector_kernel_reorder_mod.o

clean:
	rm -f *.o *.mod *.optrpt *.s kdriver_orig kdriver_xsmm kdriver_reorder

allclean: clean
