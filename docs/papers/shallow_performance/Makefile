all: figs shallow_perf.tex shallow_perf.bib
	bibtex shallow_perf
	pdflatex shallow_perf.tex

figs: best_psykal_summary.pdf orig_summary.pdf slowdown_summary.pdf

%.pdf: %.dat
	bargraph.pl $< > $*.eps
	epstopdf $*.eps
	rm $*.eps
