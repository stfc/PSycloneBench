FC=ifort
FFLAGS=-qopenmp -Ofast
#FFLAGS=-g -qopenmp -O0 -check all -traceback
LDFLAGS=
MKLFLAGS= -L${MKLROOT}/lib/intel64 -mkl
LIBXSMM=${HOME}/tools/libxsmm

VERSIONS=original fuseinterleaved kinner nlayersf nlayersf2 nlayersf_split vlen specialized specialized_save #dgemv dgemm libxsmm 

all : kdriver.avx2 kdriver.avx512 kdriver.novec
#kdriver.novec

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

matrix_vector_kernel_mod.o: matrix_vector_kernel_mod.F90
	$(FC) $(FFLAGS) -qopt_report=5 -Fa -c $<


%.novec.o: %.F90
	$(FC) -Ofast -no-vec -no-simd -qopt_report=5 -c $< -o $@

clean:
	rm -f *.o *.mod kdriver.avx512.* kdriver.novec.* *.optrpt *.s kdriver.avx2.*

tar:
	tar -cvf dino.saur Makefile *.F90
	gzip dino.saur

# dependencies
dino_mod.o: constants_mod.o
matrix_vector_kernel_mod.o: constants_mod.o argument_mod.o kernel_mod.o
kdriver.o: utils.o constants_mod.o dino_mod.o matrix_vector_kernel_mod.o

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

		#$(FC) $(FFLAGS) -no-vec -no-simd -I${HOME}/tools/libxsmm/include -DVERSION=$$v -DNEW_DL=1 -o $@.$$v kdriver.F90 $^ $(LDFLAGS) -mkl; \

run:
	@for v in $(VERSIONS); do \
		echo "Running $$v"; \
		if [ $$v = "libxsmm" ]; then \
		LIBXSMM_VERBOSE=1 ./kdriver.avx512.$$v; \
		else \
		./kdriver.avx512.$$v; \
		fi \
	done

versions:
	@echo $(VERSIONS)


