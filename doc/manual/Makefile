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

SRC = $(wildcard src/source/*.rst)

.PHONY: all html pdf

all: html pdf

html: src/build/html/index.html

src/build/html/index.html: $(SRC) src/source/conf.py
	make -C src html
	rm -rf src/build/html/_static/mathjax
	ln -s /usr/share/javascript/mathjax src/build/html/_static/mathjax

pdf: src/build/latex/dynare.pdf

src/build/latex/dynare.pdf: $(SRC) src/source/conf.py
	make -C src latexpdf

push: html
	rsync -avz src/build/html/* $(TARGET)
