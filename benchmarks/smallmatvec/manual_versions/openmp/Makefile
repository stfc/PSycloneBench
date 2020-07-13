FC=ifort
FFLAGS=-qopenmp -Ofast
#FFLAGS=-g -qopenmp -O0 -check all -traceback
LDFLAGS=

VERSIONS=original fuseinterleaved kinner nlayersf nlayersf2 nlayersf_split \
		 vlen specialized specialized_save #dgemv dgemm

# For the 'dgemv' and 'dgemm' versions the following parameters need to be
# uncommented and pointing to the appropriate library installations.
# MKLFLAGS="-L${MKLROOT}/lib/intel64 -mkl"
# LDFLAGS += ${MKLFLAGS}

all : kdriver.avx2 kdriver.avx512 kdriver.novec

# Dependencies
dino_mod.o: constants_mod.o
matrix_vector_kernel_mod.o: constants_mod.o argument_mod.o kernel_mod.o
kdriver.o: utils.o constants_mod.o dino_mod.o matrix_vector_kernel_mod.o

# Compile object targets
%.o: %.F90
	$(FC) $(FFLAGS) -c $<

matrix_vector_kernel_mod.o: matrix_vector_kernel_mod.F90
	$(FC) $(FFLAGS) -qopt_report=5 -Fa -c $<

%.novec.o: %.F90
	$(FC) $(FFLAGS) -no-vec -no-simd -qopt_report=5 -c $< -o $@

# Compile binary targets
kdriver.avx512: utils.o constants_mod.o argument_mod.o dino_mod.o kernel_mod.o matrix_vector_kernel_mod.o 
	for v in $(VERSIONS); do \
		$(FC) $(FFLAGS) -xMIC-AVX512 -DVERSION=$$v -o $@.$$v kdriver.F90 $^ $(LDFLAGS); \
	done

kdriver.avx2: utils.o constants_mod.o argument_mod.o dino_mod.o kernel_mod.o matrix_vector_kernel_mod.o 
	for v in $(VERSIONS); do \
		$(FC) $(FFLAGS) -xAVX2 -DVERSION=$$v -o $@.$$v kdriver.F90 $^ $(LDFLAGS); \
	done

kdriver.novec: utils.o constants_mod.o argument_mod.o dino_mod.o kernel_mod.o matrix_vector_kernel_mod.novec.o 
	for v in $(VERSIONS); do \
		$(FC) $(FFLAGS) -no-vec -no-simd -DVERSION=$$v -o $@.$$v kdriver.F90 $^ $(LDFLAGS) ; \
	done

run:
	@for v in $(VERSIONS); do \
		echo "Running $$v"; \
		./kdriver.avx512.$$v; \
	done

versions:
	@echo $(VERSIONS)

clean:
	rm -f *.o *.mod kdriver.avx512.* kdriver.novec.* *.optrpt *.s kdriver.avx2.*

allclean: clean

tar:
	tar -cvf dino.saur Makefile *.F90
	gzip dino.saur
