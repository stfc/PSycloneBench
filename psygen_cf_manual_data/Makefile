domain_sizes := 32 64 128 256 512 1024

%.eps: %.agr
	gracebat -hdevice EPS -printfile $@ $<

plots: $(patsubst %,sep_para_dos_%.eps, $(domain_sizes)) \
       $(patsubst %,single_region_%.eps, $(domain_sizes))

clean:
	rm *.eps

