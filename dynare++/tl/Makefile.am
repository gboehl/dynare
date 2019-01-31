SUBDIRS = cc testing

EXTRA_DIST = tl.tex

if HAVE_PDFLATEX
pdf-local: tl.pdf

tl.pdf: tl.tex
	$(PDFLATEX) tl
	$(PDFLATEX) tl
endif

CLEANFILES = *.pdf *.log *.aux
