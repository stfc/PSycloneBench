all: figs nemolite2d_perf.tex bibliography
	pdflatex nemolite2d_perf.tex

bibliography: nemolite2d_perf.bbl

nemolite2d_perf.bbl: nemolite2d_perf.aux nemolite2d_perf.bib
	bibtex nemolite2d_perf

nemolite2d_perf.aux: nemolite2d_perf.tex
	pdflatex $<

# List of figures that are GENERATED from .dat files of the same name
FIG_FILES = orig_summary.pdf \
opt_stages_256.pdf \
best_psykal_summary.pdf \
slowdown_summary.pdf \
gpu_opt_stages.pdf \
cpu_cf_gpu.pdf

figs: ${FIG_FILES}

%.pdf: %.dat
	./bargraph.pl $< > $*.eps
	epstopdf $*.eps
	rm $*.eps

clean:
	rm -f *.blg *.log *~

allclean: clean
	rm -f ${FIG_FILES}
	rm -f nemolite2d_perf.bbl nemolite2d_perf.aux nemolite2d_perf.pdf
