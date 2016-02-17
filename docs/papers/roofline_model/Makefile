all: figs gocean_roofline.tex  bibliography
	pdflatex gocean_roofline.tex

bibliography: gocean_roofline.bbl

gocean_roofline.bbl: gocean_roofline.aux gocean_roofline.bib
	bibtex gocean_roofline

gocean_roofline.aux: gocean_roofline.tex
	pdflatex $<

# List of figures that are GENERATED from .dat files of the same name
FIG_FILES = continuity_opt.pdf

# Plots generated from Grace (xmgrace) - i.e. .agr files
GRACE_PLOTS = nemolite2d_kernel_perf.pdf nemolite2d_cont_perf.pdf daxpy_bw_perf.pdf

# Plots generated from gnuplot
GNU_PLOTS = roofline_archer.pdf roofline_haswell.pdf stream2_desktop.pdf stream2_archer.pdf roofline_daxpy.pdf daxpy_plot.pdf

figs: ${FIG_FILES} ${GRACE_PLOTS} ${GNU_PLOTS}

%.pdf: %.dat
	./bargraph.pl $< > $*.eps
	epstopdf $*.eps
	rm $*.eps

%.pdf: %.agr
	gracebat -hdevice EPS -printfile $*.eps $<
	epstopdf $*.eps
	rm $*.eps

%.pdf: %.gnuplot
	gnuplot $<
	epstopdf $*.eps
	rm $*.eps

clean:
	rm -f *.blg *.log *~

allclean: clean
	rm -f ${FIG_FILES} ${GNU_PLOTS} ${GRACE_PLOTS}
	rm -f *.bbl *.aux gocean_roofline.pdf
