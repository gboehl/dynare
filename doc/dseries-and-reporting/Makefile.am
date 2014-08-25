if HAVE_PDFLATEX
if HAVE_BEAMER
pdf-local: dseriesReporting.pdf
endif
endif

SRC = dseriesReporting.tex

EXTRA_DIST = $(SRC)

dseriesReporting.pdf: $(SRC)
	$(PDFLATEX) dseriesReporting
	$(PDFLATEX) dseriesReporting

clean-local:
	rm -f dseriesReporting.pdf *.toc *.aux *.log *.nav *.snm *.vrb *.out *~
