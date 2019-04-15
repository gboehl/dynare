SUBDIRS = cc testing

EXTRA_DIST = sylvester.tex change_log.html

if HAVE_PDFLATEX
pdf-local: sylvester.pdf

sylvester.pdf: sylvester.tex
	$(PDFLATEX) sylvester
	$(PDFLATEX) sylvester
endif

CLEANFILES = *.pdf *.log *.aux

