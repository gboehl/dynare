SUBDIRS = cc testing

PDF_TARGETS =

if HAVE_PDFLATEX
PDF_TARGETS += tl.pdf
endif

pdf-local: $(PDF_TARGETS)

tl.pdf: tl.tex
	$(PDFLATEX) tl
	$(PDFLATEX) tl

clean-local:
	rm -f *~ *.pdf *.log *.aux *.out
