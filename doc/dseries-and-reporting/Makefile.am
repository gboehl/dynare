if ENABLE_DOC
pdf-local: dseriesReporting.pdf
endif

SRC = dseriesReporting.tex

EXTRA_DIST = $(SRC)

dseriesReporting.pdf: $(SRC)
	$(PDFLATEX) dseriesReporting
	$(PDFLATEX) dseriesReporting

clean-local:
	rm -f dseriesReporting.pdf *.toc *.aux *.log *.nav *.snm *.vrb *.out *~
