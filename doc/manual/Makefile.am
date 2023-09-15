EXTRA_DIST = source \
	utils/dynare_dom.py \
	utils/dynare_lex.py

SRC = $(wildcard source/*.rst)

html-local: build/html/index.html

build/html/index.html: $(SRC) source/conf.py
	$(SPHINXBUILD) -b html -D version=$(PACKAGE_VERSION) -D release=$(PACKAGE_VERSION) -d build/doctrees source build/html

pdf-local: build/latex/dynare-manual.pdf

# NB: we cannot use -M latexpdf because of this bug:
#  https://bugs.debian.org/933347
build/latex/dynare-manual.pdf: $(SRC) source/conf.py
	$(SPHINXBUILD) -b latex -D version=$(PACKAGE_VERSION) -D release=$(PACKAGE_VERSION) -d build/doctrees source build/latex
	make -C build/latex all-pdf

clean-local:
	rm -rf build
	rm -rf utils/__pycache__
