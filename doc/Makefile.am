SUBDIRS = parallel internals gsa dseries-and-reporting

if ENABLE_DOC
SUBDIRS += manual
pdf-local: guide.pdf bvar-a-la-sims.pdf dr.pdf sylvester.pdf tl.pdf
endif

EXTRA_DIST = guide.tex bvar-a-la-sims.tex dr.tex dr.bib dynare.plots

guide.pdf: guide.tex
	$(PDFLATEX) guide
	$(PDFLATEX) guide

bvar-a-la-sims.pdf: bvar-a-la-sims.tex
	$(PDFLATEX) bvar-a-la-sims
	$(PDFLATEX) bvar-a-la-sims

dr.pdf: dr.tex dr.bib
	$(PDFLATEX) dr
	$(BIBTEX) dr
	$(PDFLATEX) dr
	$(PDFLATEX) dr

sylvester.pdf: sylvester.tex
	$(PDFLATEX) sylvester
	$(PDFLATEX) sylvester

tl.pdf: tl.tex
	$(PDFLATEX) tl
	$(PDFLATEX) tl

clean-local:
	rm -f *~ *.pdf *.log *.aux *.out *.blg *.bbl
