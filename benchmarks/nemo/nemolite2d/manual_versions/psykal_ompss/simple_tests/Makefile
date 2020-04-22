F90 = gfortran
F90FLAGS =
#F90 = mfc
#F90FLAGS = --ompss
FOBJ = field.o main.o psy.o my_work1.o my_work2.o

test_example: ${FOBJ}
	${F90} -o $@ ${FOBJ} ${F90FLAGS}

psy.o : my_work1.o my_work2.o field.o
main.o : psy.o field.o

%.o: %.f90
	${F90} ${F90FLAGS} -c $<
clean:
	rm -f *.o *.mod test_example
