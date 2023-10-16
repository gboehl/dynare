# Copyright © 2019-2023 Dynare Team
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
# along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

ARCHS = arm64 x86_64
BUILDS = $(addprefix build-, $(ARCHS))
DEPS = $(addprefix deps-, $(ARCHS))
CLEAN_DEPS = $(addprefix clean-deps-, $(ARCHS))
CLEAN_ALL = $(addprefix clean-all-, $(ARCHS))

.PHONY: all $(BUILDS) $(DEPS) clean-pkg $(CLEAN_DEPS) $(CLEAN_ALL) clean-all

all: $(BUILDS)

$(BUILDS): build-%: deps-%
	./build.sh $*

$(DEPS): deps-%:
	make -C deps DEPS_ARCH=$*

$(CLEAN_DEPS): clean-deps-%:
	make -C deps clean-all DEPS_ARCH=$*

clean-pkg:
	rm -rf pkg

$(CLEAN_ALL): clean-all-%: clean-deps-% clean-pkg

clean-all: $(CLEAN_ALL)
