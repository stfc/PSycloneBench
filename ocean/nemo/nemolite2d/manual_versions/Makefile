# This Makefile iterates through the IMPLEMENTATION directories and
# compiles, executes and prints a summary with a copy of an equivalent
# namelist, or cleans, each directory.


# All compilers should compile this IMPLEMENTATIONS
IMPLEMENTATIONS = psykal_serial psykal_omp
# Uncomment IMPLEMENTATIONS below if your environment has MPI compilers
# IMPLEMENTATIONS += psykal_dm
# Uncomment IMPLEMENTATIONS below if your environment has OpenCL compilers
# IMPLEMENTATIONS += psykal_ocl
# Uncomment IMPLEMENTATIONS below if your environment has OMPSS compilers
# IMPLEMENTATIONS += psykal_ompss

all: compile

# Compile all enabled implementations
compile:
	$(foreach var,$(IMPLEMENTATIONS), \
		cd $(var) && $(MAKE) nemolite2d; cd ..;)

# Execute all enabled implementations with their best full node performance and
# report their ua and uv final checksums and the time-stepping time.
run:
	@echo "Implementation\t\t ua checksum\t uv checksum\t time-stepping"
	@$(foreach var,$(IMPLEMENTATIONS), \
		cp namelist $(var)/. && \
		cd $(var) && \
		echo -n "$(var): \t\t" && \
		$(MAKE) --no-print-directory summary && \
		cd ..; \
	)

# Clean folders of all enabled implementations
clean:
	$(foreach var,$(IMPLEMENTATIONS),cd $(var) && $(MAKE) clean && cd ..;)
