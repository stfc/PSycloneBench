# This Makefile iterates through the IMPLEMENTATION directories and
# compiles, executes and prints a summary with a copy of an equivalent
# namelist, or cleans, each directory.


IMPLEMENTATIONS = psykal_serial psykal_omp

all: compile

compile:
	$(foreach var,$(IMPLEMENTATIONS),cd $(var) && $(MAKE) nemolite2d && cd ..;)

run:
	@echo "Implementation\t ua checksum\t uv checksum\t time-stepping"
	@$(foreach var,$(IMPLEMENTATIONS), \
		cp namelist $(var)/. && \
		cd $(var) && \
		echo -n "$(var): " && \
		$(MAKE) --no-print-directory summary && \
		cd ..; \
	)

clean:
	$(foreach var,$(IMPLEMENTATIONS),cd $(var) && $(MAKE) clean && cd ..;)
