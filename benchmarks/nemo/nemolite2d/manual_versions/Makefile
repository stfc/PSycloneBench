# This Makefile iterates through the IMPLEMENTATION directories and
# either cleans or compiles, executes and prints a summary with a copy
# of an equivalent namelist in each one.

# All compilers should be able to compile these IMPLEMENTATIONS
IMPLEMENTATIONS = psykal_serial
IMPLEMENTATIONS += psykal_cpp
IMPLEMENTATIONS += psykal_omp
# Uncomment IMPLEMENTATIONS below if your environment has MPI compilers
# IMPLEMENTATIONS += psykal_dm
# Uncomment IMPLEMENTATIONS below if your environment has the Kokkos libraries
IMPLEMENTATIONS += psykal_kokkos
# Uncomment IMPLEMENTATIONS below if your environment has OpenCL compilers
IMPLEMENTATIONS += psykal_opencl psykal_c_opencl
# Uncomment IMPLEMENTATIONS below if your environment has the mcxx OmpSs compiler.
# IMPLEMENTATIONS += psykal_cpp_ompss
# psykal_ompss and psykal_acc have not been adapted to this Makefile format
# yet.

.PHONY: all clean allclean compile summary_table
all: compile

# Compile all enabled implementations (some targets need to start with an
# allclean)
compile:
	$(foreach var, $(IMPLEMENTATIONS), \
		cd $(var) && $(MAKE) allclean && $(MAKE); cd ..;)

# Execute all enabled implementations with their best full node performance and
# report their ua and uv final checksums and the time-stepping time.
summary_table:
	@echo "Implementation\t\tCompiler\tArch\tua checksum\tuv checksum\ttime/step"
	@$(foreach var, $(IMPLEMENTATIONS), \
		cd $(var) && \
		echo -n "$(var): \t\t" && \
		$(MAKE) --no-print-directory summary && \
		cd ..; \
	)

# Clean folders (and dependencies) of all enabled implementations
clean:
	$(foreach var, $(IMPLEMENTATIONS), cd $(var) && $(MAKE) clean && cd ..;)

allclean:
	$(foreach var, $(IMPLEMENTATIONS), cd $(var) && $(MAKE) allclean && cd ..;)
