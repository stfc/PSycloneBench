

psy_layer: psy_layer.F90
	${F90} ${F90FLAGS} -c $<

openmp_cpu: psy_layer.F90
	${F90} ${F90FLAGS} ${OMPFLAGS} -c $< -o openmp_cpu_psy_layer.o
	
openmp_gpu: psy_layer.F90
	${F90} ${F90FLAGS} -DTARGET_GPU ${OMPTARGETFLAGS} ${UMEMFLAGS} -c $< -o openmp_gpu_psy_layer.o

clean:
	rm -rf *.o *.mod
