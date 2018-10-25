all: html pdf

html: src/build/html/index.html

src/build/html/index.html: $(src/source/%.rst) src/source/conf.py
	source python/bin/activate ; make -C src html

pdf: src/build/latex/dynare.pdf

src/build/latex/dynare.pdf: $(src/source/%.rst) src/source/conf.py
	source python/bin/activate ; make -C src latexpdf

python: python/bin/python3

python/bin/python3:
	python3 -m venv python
	source python/bin/activate ; pip install --upgrade pip ; pip install sphinx recommonmark sphinx_rtd_theme
	cp py/pygment/dynare.py python/lib/python3.*/site-packages/pygments/lexers/
	cd python/lib/python3.*/site-packages/pygments/lexers ; python3 _mapping.py
