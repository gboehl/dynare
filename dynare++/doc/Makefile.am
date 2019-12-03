EXTRA_DIST = \
	dynare++-ramsey.tex \
	dynare++-tutorial.tex \
	sylvester.tex \
	tl.tex \
	changelog-old.html \
	changelog-sylv-old.html

if ENABLE_DOC
pdf-local: dynare++-ramsey.pdf dynare++-tutorial.pdf sylvester.pdf tl.pdf
endif

%.pdf: %.tex
	$(PDFLATEX) $<
	$(PDFLATEX) $<
	$(PDFLATEX) $<

CLEANFILES = *.pdf *.log *.aux *.out *.toc


