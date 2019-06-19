SUBDIRS = cc testing

EXTRA_DIST = sylvester.tex

if HAVE_PDFLATEX
pdf-local: sylvester.pdf

sylvester.pdf: sylvester.tex
	$(PDFLATEX) sylvester
	$(PDFLATEX) sylvester
endif

CLEANFILES = *.pdf *.log *.aux

