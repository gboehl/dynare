# Copyright (C) 2018-2019 Dynare Team
#
# This file is part of Dynare.
#
# Dynare is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Dynare is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

MATHJAX_VERSION = 2.7.5

SRC = $(wildcard src/source/*.rst)

.PHONY: all html pdf deps python mathjax

all: html pdf

html: deps src/build/html/index.html

src/build/html/index.html: $(SRC) src/source/conf.py
	. python/bin/activate ; make -C src html

pdf: deps src/build/latex/dynare.pdf

src/build/latex/dynare.pdf: $(SRC) src/source/conf.py
	. python/bin/activate ; make -C src latexpdf

deps: python mathjax

python: python/bin/python3

python/bin/python3:
	python3 -m venv python
	. python/bin/activate ; pip3 install --upgrade pip ; pip3 install sphinx recommonmark sphinx_rtd_theme
	cp py/pygment/dynare.py python/lib/python3.*/site-packages/pygments/lexers/
	cd python/lib/python3.*/site-packages/pygments/lexers ; python3 _mapping.py

mathjax: src/source/_static/mathjax/MathJax.js
	@touch src/source/_static/mathjax/MathJax.js

src/source/_static/mathjax/MathJax.js: mathjax-$(MATHJAX_VERSION).zip
	unzip mathjax-$(MATHJAX_VERSION).zip
	mv MathJax-$(MATHJAX_VERSION) src/source/_static/mathjax

mathjax-$(MATHJAX_VERSION).zip:
	wget https://github.com/mathjax/MathJax/archive/$(MATHJAX_VERSION).zip
	mv $(MATHJAX_VERSION).zip mathjax-$(MATHJAX_VERSION).zip
	@touch mathjax-$(MATHJAX_VERSION).zip
