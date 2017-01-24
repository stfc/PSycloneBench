all: figs habakkuk.tex  bibliography
	pdflatex habakkuk.tex

bibliography: habakkuk.bbl

habakkuk.bbl: habakkuk.aux habakkuk.bib
	bibtex habakkuk

habakkuk.aux: habakkuk.tex
	pdflatex $<

# List of figures that are GENERATED from .dat files of the same name
FIG_FILES = mulsd_divsd.pdf

# Plots generated from Grace (xmgrace) - i.e. .agr files
GRACE_PLOTS = 

# Plots generated from gnuplot
GNU_PLOTS =

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
