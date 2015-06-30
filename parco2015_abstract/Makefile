BASE_NAME = psyclone_parco2015

all: figs ${BASE_NAME}.tex bibliography
	pdflatex ${BASE_NAME}.tex

bibliography: ${BASE_NAME}.bbl

${BASE_NAME}.bbl: ${BASE_NAME}.aux ../shallow_perf.bib
	bibtex ${BASE_NAME}

${BASE_NAME}.aux: ${BASE_NAME}.tex
	pdflatex $<

# List of figures that are GENERATED from .dat files of the same name
FIG_FILES = slowdown_summary.pdf

figs: ${FIG_FILES}

%.pdf: %.dat
	../bargraph.pl $< > $*.eps
	epstopdf $*.eps
	rm $*.eps

clean:
	rm -f *.blg *.log *~

allclean: clean
	rm -f ${FIG_FILES}
	rm -f ${BASE_NAME}.bbl ${BASE_NAME}.aux ${BASE_NAME}.pdf
