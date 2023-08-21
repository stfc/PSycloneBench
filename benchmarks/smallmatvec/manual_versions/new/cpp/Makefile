

psy_layer: psy_layer.f90 cpp_psy_layer.o
	${F90} ${FFLAGS} -c $<

cpp_psy_layer.o: psy_layer.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@


clean:
	rm -rf *.o *.mod
