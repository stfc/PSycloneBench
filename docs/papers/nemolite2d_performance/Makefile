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

# Plots generated from Grace (xmgrace)
GRACE_PLOTS = omp_scaling_32_gnu.pdf \
              omp_scaling_32_cray.pdf \
              omp_scaling_32_intel.pdf \
              omp_scaling_problem_size.pdf

figs: ${FIG_FILES} ${GRACE_PLOTS}

%.pdf: %.dat
	./bargraph.pl $< > $*.eps
	epstopdf $*.eps
	rm $*.eps

%.pdf: %.agr
	gracebat -hdevice EPS -printfile $*.eps $<
	epstopdf $*.eps
	rm $*.eps

clean:
	rm -f *.blg *.log *~

allclean: clean
	rm -f ${FIG_FILES}
	rm -f nemolite2d_perf.bbl nemolite2d_perf.aux nemolite2d_perf.pdf
