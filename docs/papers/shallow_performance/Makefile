all: figs shallow_perf.tex bibliography
	pdflatex shallow_perf.tex

bibliography: shallow_perf.bbl

shallow_perf.bbl: shallow_perf.aux shallow_perf.bib
	bibtex shallow_perf

shallow_perf.aux: shallow_perf.tex
	pdflatex $<

# List of figures that are GENERATED from .dat files of the same name
FIG_FILES = best_psykal_summary.pdf orig_summary.pdf slowdown_summary.pdf opt_stages_256.pdf

figs: ${FIG_FILES}

%.pdf: %.dat
	./bargraph.pl $< > $*.eps
	epstopdf $*.eps
	rm $*.eps

clean:
	rm -f *.blg *.log *~

allclean: clean
	rm -f ${FIG_FILES}
	rm -f shallow_perf.bbl shallow_perf.aux shallow_perf.pdf
