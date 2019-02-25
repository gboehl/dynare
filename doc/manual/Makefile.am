EXTRA_DIST = source \
	utils/dynare_dom.py \
	utils/dynare_lex.py

SRC = $(wildcard src/source/*.rst)

html-local: build/html/index.html

build/html/index.html: $(SRC) source/conf.py
	$(SPHINXBUILD) -M html source build

pdf-local: build/latex/dynare-manual.pdf

build/latex/dynare-manual.pdf: $(SRC) source/conf.py
	$(SPHINXBUILD) -M latexpdf source build

clean-local:
	rm -rf build
