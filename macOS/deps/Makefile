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

include versions.mk

# settings for different architectures
DEPS_ARCH ?= x86_64 # use x86_64 by default

BREWDIR := $(if $(filter arm64,$(DEPS_ARCH)),/opt/homebrew,/usr/local)

GCC_VERSION = 13
ROOT_PATH = $(realpath .)

WGET_OPTIONS := --no-verbose --no-use-server-timestamps --retry-connrefused --retry-on-host-error

NTHREADS=$(shell sysctl -n hw.perflevel0.physicalcpu)

.PHONY: all build build-slicot build-x13as \
	clean-lib clean-libslicot clean-x13as-bin \
	clean-src clean-slicot-src clean-x13as-src \
	clean-tar clean-slicot-tar clean-x13as-tar \
	clean-all clean-lib clean-src clean-tar

all: build

build: build-slicot build-x13as ln-matio ln-gsl

clean-lib: clean-libslicot clean-x13as-bin

clean-src: clean-slicot-src clean-x13as-src

clean-tar: clean-slicot-tar clean-x13as-tar

clean-all: clean-lib clean-src clean-tar


#
# Matio & GSL
# (done to link only to static Matio and GSL libraries)

# Matio
$(DEPS_ARCH)/lib64/matio/lib/libmatio.a: $(BREWDIR)/lib/libmatio.a
	mkdir -p $(dir $@) && ln -sf $< $@

$(DEPS_ARCH)/lib64/matio/lib/libhdf5.a: $(BREWDIR)/lib/libhdf5.a
	mkdir -p $(dir $@) && ln -sf $< $@

$(DEPS_ARCH)/lib64/matio/lib/libsz.a: $(BREWDIR)/lib/libsz.a
	mkdir -p $(dir $@) && ln -sf $< $@

$(DEPS_ARCH)/lib64/matio/include/matio.h: $(BREWDIR)/include/matio.h
	mkdir -p $(dir $@) && cd $(dir $@).. && rm -rf include && ln -sf $(dir $<) .

ln-matio: $(DEPS_ARCH)/lib64/matio/lib/libmatio.a \
	$(DEPS_ARCH)/lib64/matio/lib/libhdf5.a \
	$(DEPS_ARCH)/lib64/matio/lib/libsz.a \
	$(DEPS_ARCH)/lib64/matio/include/matio.h

clean-matio:
	rm -rf $(DEPS_ARCH)/lib64/matio

# GSL
$(DEPS_ARCH)/lib64/gsl/lib/libgsl.a: $(BREWDIR)/lib/libgsl.a
	mkdir -p $(dir $@) && ln -sf $< $@

$(DEPS_ARCH)/lib64/gsl/lib/libgslcblas.a: $(BREWDIR)/lib/libgslcblas.a
	mkdir -p $(dir $@) && ln -sf $< $@

$(DEPS_ARCH)/lib64/gsl/include/gsl/gsl_blas.h: $(BREWDIR)/include/gsl/gsl_blas.h
	mkdir -p $(dir $@) && cd $(dir $@).. && rm -rf gsl && ln -sf $(dir $<) .

ln-gsl: $(DEPS_ARCH)/lib64/gsl/lib/libgsl.a \
	$(DEPS_ARCH)/lib64/gsl/lib/libgslcblas.a \
	$(DEPS_ARCH)/lib64/gsl/include/gsl/gsl_blas.h

clean-gsl:
	rm -rf $(DEPS_ARCH)/lib64/gsl

#
# slicot
#
tarballs/slicot-$(SLICOT_VERSION).tar.gz:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://deb.debian.org/debian/pool/main/s/slicot/slicot_$(SLICOT_VERSION).orig.tar.gz

$(DEPS_ARCH)/sources64/slicot-$(SLICOT_VERSION)-with-32bit-integer-and-underscore: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	rm -rf $(DEPS_ARCH)/sources64/slicot-*-with-32bit-integer-and-underscore
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1
	touch $@

$(DEPS_ARCH)/sources64/slicot-$(SLICOT_VERSION)-with-64bit-integer-and-underscore: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	rm -rf $(DEPS_ARCH)/sources64/slicot-*-with-64bit-integer-and-underscore
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1
	touch $@

$(DEPS_ARCH)/lib64/slicot/with-underscore/lib/libslicot_pic.a: $(DEPS_ARCH)/sources64/slicot-$(SLICOT_VERSION)-with-32bit-integer-and-underscore
	make -C $< FORTRAN=$(BREWDIR)/bin/gfortran LOADER=$(BREWDIR)/bin/gfortran SLICOTLIB=../libslicot_pic.a OPTS="-O2 -g" lib -j$(NTHREADS)
	strip -S $</libslicot_pic.a
	mkdir -p $(dir $@)
	cp $</libslicot_pic.a $@

$(DEPS_ARCH)/lib64/slicot/with-underscore/lib/libslicot64_pic.a: $(DEPS_ARCH)/sources64/slicot-$(SLICOT_VERSION)-with-64bit-integer-and-underscore
	make -C $< FORTRAN=$(BREWDIR)/bin/gfortran LOADER=$(BREWDIR)/bin/gfortran SLICOTLIB=../libslicot64_pic.a OPTS="-O2 -g -fdefault-integer-8" lib -j$(NTHREADS)
	strip -S $</libslicot64_pic.a
	mkdir -p $(dir $@)
	cp $</libslicot64_pic.a $@

build-slicot: $(DEPS_ARCH)/lib64/slicot/with-underscore/lib/libslicot_pic.a \
	$(DEPS_ARCH)/lib64/slicot/with-underscore/lib/libslicot64_pic.a

clean-slicot-tar:
	rm -f tarballs/slicot-$(SLICOT_VERSION).tar.gz

clean-slicot-src:
	rm -rf $(DEPS_ARCH)/sources64/slicot-$(SLICOT_VERSION)-with-64bit-integer-and-underscore

clean-libslicot:
	rm -rf $(DEPS_ARCH)/lib64/Slicot

clean-slicot-all: clean-slicot-src clean-slicot-tar clean-libslicot



#
# X13AS
#
tarballs/x13as_asciisrc-v$(X13AS_VERSION).tar.gz:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://www2.census.gov/software/x-13arima-seats/x13as/unix-linux/program-archives/x13as_asciisrc-v$(X13AS_VERSION).tar.gz

$(DEPS_ARCH)/sources64/x13as-$(X13AS_VERSION): tarballs/x13as_asciisrc-v$(X13AS_VERSION).tar.gz
	rm -rf $(DEPS_ARCH)/sources64/x13as-*
	mkdir -p $@
	tar xf $< --directory $@

$(DEPS_ARCH)/lib64/x13as/x13as: $(DEPS_ARCH)/sources64/x13as-$(X13AS_VERSION)
	# Statically link x13as (see #1865).
	# Using -static is not possible, it does not work under Darwin.
	# Implement an ugly workaround in the absence of -static-libquadmath flag.
	# In particular, gcc must be used as the linker, because gfortran does not
	# honour static linking of libquadmath.
	# Once the -static-libquadmath flag is added (GCC 13? ; see
	# https://gcc.gnu.org/bugzilla/show_bug.cgi?id=46539), go back to using
	# gfortran as the linker with -static-libgfortran and
	# -static-libquadmath flags, and drop the GCC_VERSION variable.
	cd $< && sed -i '' 's/-static//g' makefile.gf
	make -C $< -f makefile.gf FC=$(BREWDIR)/bin/gfortran LINKER=$(BREWDIR)/bin/gcc-$(GCC_VERSION) FFLAGS="-O2 -std=legacy" LDFLAGS=-static-libgcc LIBS="$(BREWDIR)/lib/gcc/current/libgfortran.a $(BREWDIR)/lib/gcc/current/libquadmath.a" PROGRAM=x13as -j$(NTHREADS)
	strip $</x13as
	mkdir -p $(dir $@)
	cp $</x13as $@

build-x13as: $(DEPS_ARCH)/lib64/x13as/x13as

clean-x13as-tar:
	rm -f tarballs/x13as_asciisrc-v$(X13AS_VERSION).tar.gz

clean-x13as-src:
	rm -rf $(DEPS_ARCH)/sources64/x13as-$(X13AS_VERSION)

clean-x13as-bin:
	rm -rf $(DEPS_ARCH)/lib64/x13as

clean-x13as-all: clean-x13as-tar clean-x13as-src clean-x13as-bin
