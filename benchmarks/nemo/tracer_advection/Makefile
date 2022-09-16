.PHONY: all compute_in_subroutine multi_kernel original

all: compute_in_subroutine multi_kernel original

compute_in_subroutine multi_kernel original:
	${MAKE} -C $@

clean allclean:
	${MAKE} -C compute_in_subroutine $@
	${MAKE} -C multi_kernel $@
	${MAKE} -C original $@
