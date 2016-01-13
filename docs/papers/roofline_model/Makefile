all: figs gocean_roofline.tex  bibliography
	pdflatex gocean_roofline.tex

bibliography: gocean_roofline.bbl

gocean_roofline.bbl: gocean_roofline.aux gocean_roofline.bib
	bibtex gocean_roofline

gocean_roofline.aux: gocean_roofline.tex
	pdflatex $<

# List of figures that are GENERATED from .dat files of the same name
FIG_FILES = 

# Plots generated from Grace (xmgrace)
GRACE_PLOTS = 

# Plots generated from gnuplot
GNU_PLOTS = roofline_archer.pdf roofline_haswell.pdf

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
	rm -f ${FIG_FILES}
	rm -f *.bbl *.aux gocean_roofline.pdf
