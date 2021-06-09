# Copyright © 2017-2020 Dynare Team
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

ROOT_PATH = $(realpath .)

WGET_OPTIONS := --no-verbose --no-use-server-timestamps --retry-connrefused --retry-on-host-error

.PHONY: all build octave matlab msys2 \
	build build-slicot build-x13as \
	clean-lib clean-libslicot clean-matlab clean-octave clean-msys2 clean-x13as-bin \
	clean-src clean-slicot-src clean-x13as-src \
	clean-tar clean-slicot-tar clean-x13as-tar \
	clean-all clean-lib clean-src clean-tar

all: build octave matlab msys2

build: build-slicot build-x13as

clean-lib: clean-libslicot clean-matlab clean-octave clean-msys2 clean-x13as-bin

clean-src: clean-slicot-src clean-x13as-src

clean-tar: clean-slicot-tar clean-x13as-tar

clean-all: clean-lib clean-src clean-tar

#
# Slicot
#

tarballs/slicot-$(SLICOT_VERSION).tar.gz:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://deb.debian.org/debian/pool/main/s/slicot/slicot_$(SLICOT_VERSION).orig.tar.gz

sources64/slicot-$(SLICOT_VERSION)-with-32bit-integer: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	rm -rf sources64/slicot-*-with-32bit-integer
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1
	touch $@

sources64/slicot-$(SLICOT_VERSION)-with-64bit-integer: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	rm -rf sources64/slicot-*-with-64bit-integer
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1
	touch $@

sources64/slicot-$(SLICOT_VERSION)-with-32bit-integer-and-underscore: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	rm -rf sources64/slicot-*-with-32bit-integer-and-underscore
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1
	touch $@

sources64/slicot-$(SLICOT_VERSION)-with-64bit-integer-and-underscore: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	rm -rf sources64/slicot-*-with-64bit-integer-and-underscore
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1
	touch $@

lib64/Slicot/without-underscore/lib/libslicot_pic.a: sources64/slicot-$(SLICOT_VERSION)-with-32bit-integer
	make -C $< lib SLICOTLIB=../libslicot_pic.a OPTS="-O2 -g -fno-underscoring" FORTRAN=x86_64-w64-mingw32-gfortran LOADER=x86_64-w64-mingw32-gfortran ARCH=x86_64-w64-mingw32-ar
	x86_64-w64-mingw32-strip --strip-debug $</libslicot_pic.a
	mkdir -p $(dir $@)
	cp $</libslicot_pic.a $@

lib64/Slicot/without-underscore/lib/libslicot64_pic.a: sources64/slicot-$(SLICOT_VERSION)-with-64bit-integer
	make -C $< lib SLICOTLIB=../libslicot64_pic.a OPTS="-O2 -g -fno-underscoring -fdefault-integer-8" FORTRAN=x86_64-w64-mingw32-gfortran LOADER=x86_64-w64-mingw32-gfortran ARCH=x86_64-w64-mingw32-ar
	x86_64-w64-mingw32-strip --strip-debug $</libslicot64_pic.a
	mkdir -p $(dir $@)
	cp $</libslicot64_pic.a $@

lib64/Slicot/with-underscore/lib/libslicot_pic.a: sources64/slicot-$(SLICOT_VERSION)-with-32bit-integer-and-underscore
	make -C $< lib SLICOTLIB=../libslicot_pic.a OPTS="-O2 -g" FORTRAN=x86_64-w64-mingw32-gfortran LOADER=x86_64-w64-mingw32-gfortran ARCH=x86_64-w64-mingw32-ar
	x86_64-w64-mingw32-strip --strip-debug $</libslicot_pic.a
	mkdir -p $(dir $@)
	cp $</libslicot_pic.a $@

lib64/Slicot/with-underscore/lib/libslicot64_pic.a: sources64/slicot-$(SLICOT_VERSION)-with-64bit-integer-and-underscore
	make -C $< lib SLICOTLIB=../libslicot64_pic.a OPTS="-O2 -g -fdefault-integer-8" FORTRAN=x86_64-w64-mingw32-gfortran LOADER=x86_64-w64-mingw32-gfortran ARCH=x86_64-w64-mingw32-ar
	x86_64-w64-mingw32-strip --strip-debug $</libslicot64_pic.a
	mkdir -p $(dir $@)
	cp $</libslicot64_pic.a $@

build-slicot: lib64/Slicot/without-underscore/lib/libslicot_pic.a \
	lib64/Slicot/without-underscore/lib/libslicot64_pic.a \
	lib64/Slicot/with-underscore/lib/libslicot_pic.a \
	lib64/Slicot/with-underscore/lib/libslicot64_pic.a

clean-slicot-tar:
	rm -f tarballs/slicot-$(SLICOT_VERSION).tar.gz

clean-slicot-src:
	rm -rf sources64/slicot-$(SLICOT_VERSION)-with-32bit-integer
	rm -rf sources64/slicot-$(SLICOT_VERSION)-with-64bit-integer
	rm -rf sources64/slicot-$(SLICOT_VERSION)-with-32bit-integer-and-underscore
	rm -rf sources64/slicot-$(SLICOT_VERSION)-with-64bit-integer-and-underscore

clean-libslicot:
	rm -rf lib64/Slicot

clean-slicot-all: clean-slicot-src clean-slicot-tar clean-libslicot

#
# Octave
#

tarballs/octave-${OCTAVE_VERSION}${OCTAVE_W64_BUILD}-w64.7z:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://ftp.gnu.org/gnu/octave/windows/$(notdir $@)

octave64: tarballs/octave-${OCTAVE_VERSION}${OCTAVE_W64_BUILD}-w64.7z
	rm -rf $@
	7zr x $< > /dev/null
	mv octave-${OCTAVE_VERSION}-w64 $@
	# Clean up MinGW stuff that interferes with our cross-compiler
	cd $@/mingw64/lib && rm -f libadvapi32.a libkernel32.a libmingw32.a libmingwex.a libmoldname.a libmsvcrt.a libpthread.dll.a libshell32.a libuser32.a
	cd $@/mingw64/include && rm -f $$(grep -l "This file is part of the mingw-w64 runtime package" *.h)
	touch $@

octave: octave64

clean-octave:
	rm -rf octave64
	rm -f tarballs/octave-${OCTAVE_VERSION}${OCTAVE_W64_BUILD}-w64.7z

#
# MATLAB
#

tarballs/matlab%.tar.xz:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@.gpg https://www.dynare.org/windows-pkg-build/$(notdir $@).gpg
	gpg --output $@ --decrypt $@.gpg
	rm $@.gpg

matlab64: tarballs/matlab64-$(MATLAB64_VERSION).tar.xz
	rm -rf $@
	mkdir -p $@
	tar xf $< --directory $@
	touch $@

matlab: matlab64

clean-matlab:
	rm -f tarballs/matlab64.tar.xz
	rm -rf matlab64

#
# Stuff downloaded from MSYS2:
# - several build dependencies
# - the embedded MinGW compiler (for use_dll)
#

msys2: lib64-msys2 mingw64

MSYS2_EXCLUDES := --exclude .MTREE --exclude .BUILDINFO --exclude .PKGINFO

lib64-msys2: tarballs/mingw-w64-x86_64-boost-$(MINGW64_BOOST_VERSION)-any.pkg.tar.zst tarballs/mingw-w64-x86_64-gsl-$(MINGW64_GSL_VERSION)-any.pkg.tar.xz tarballs/mingw-w64-x86_64-openblas-$(MINGW64_OPENBLAS_VERSION)-any.pkg.tar.zst tarballs/mingw-w64-x86_64-matio-$(MINGW64_MATIO_VERSION)-any.pkg.tar.zst tarballs/mingw-w64-x86_64-zlib-$(MINGW64_ZLIB_VERSION)-any.pkg.tar.zst tarballs/mingw-w64-x86_64-hdf5-$(MINGW64_HDF5_VERSION)-any.pkg.tar.zst tarballs/mingw-w64-x86_64-szip-$(MINGW64_SZIP_VERSION)-any.pkg.tar.xz
	rm -rf $@
	mkdir $@
	for f in $^; do tar xf $$f --directory $@ --strip-components 1 $(MSYS2_EXCLUDES); done
	touch $@

mingw64: tarballs/mingw-w64-x86_64-gcc-$(MINGW64_GCC_VERSION)-any.pkg.tar.zst tarballs/mingw-w64-x86_64-gcc-libs-$(MINGW64_GCC_VERSION)-any.pkg.tar.zst tarballs/mingw-w64-x86_64-gmp-$(MINGW64_GMP_VERSION)-any.pkg.tar.zst tarballs/mingw-w64-x86_64-binutils-$(MINGW64_BINUTILS_VERSION)-any.pkg.tar.zst tarballs/mingw-w64-x86_64-headers-git-$(MINGW64_HEADERS_VERSION)-any.pkg.tar.zst tarballs/mingw-w64-x86_64-crt-git-$(MINGW64_CRT_VERSION)-any.pkg.tar.zst tarballs/mingw-w64-x86_64-winpthreads-git-$(MINGW64_WINPTHREADS_VERSION)-any.pkg.tar.zst tarballs/mingw-w64-x86_64-libwinpthread-git-$(MINGW64_WINPTHREADS_VERSION)-any.pkg.tar.zst tarballs/mingw-w64-x86_64-zlib-$(MINGW64_ZLIB_VERSION)-any.pkg.tar.zst tarballs/mingw-w64-x86_64-zstd-$(MINGW64_ZSTD_VERSION)-any.pkg.tar.zst
	rm -rf $@
	for f in $^; do tar xf $$f $(MSYS2_EXCLUDES); done
	touch $@

tarballs/mingw-w64-x86_64-%-any.pkg.tar.xz:
	wget $(WGET_OPTIONS) -O $@ http://repo.msys2.org/mingw/x86_64/$(notdir $@)

tarballs/mingw-w64-x86_64-%-any.pkg.tar.zst:
	wget $(WGET_OPTIONS) -O $@ http://repo.msys2.org/mingw/x86_64/$(notdir $@)

clean-msys2:
	rm -rf lib64-msys2
	rm -rf mingw64
	rm -f tarballs/mingw-w64-x86_64-*.tar.xz
	rm -f tarballs/mingw-w64-x86_64-*.tar.zst

#
# X13AS
#

tarballs/x13assrc_V$(X13AS_VERSION).tar.gz:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://www.census.gov/ts/x13as/unix/x13assrc_V$(X13AS_VERSION).tar.gz

sources64/x13as-$(X13AS_VERSION): tarballs/x13assrc_V$(X13AS_VERSION).tar.gz
	rm -rf sources64/x13as-*
	mkdir -p $@
	tar xf $< --directory $@

lib64/x13as/x13as.exe: sources64/x13as-$(X13AS_VERSION)
	make -C $< -f makefile.gf FC=x86_64-w64-mingw32-gfortran LINKER=x86_64-w64-mingw32-gfortran FFLAGS="-O2 -std=legacy" PROGRAM=x13as.exe
	x86_64-w64-mingw32-strip $</x13as.exe
	mkdir -p $(dir $@)
	cp $</x13as.exe $@

build-x13as: lib64/x13as/x13as.exe

clean-x13as-tar:
	rm -f tarballs/x13assrc_V$(X13AS_VERSION).tar.gz

clean-x13as-src:
	rm -rf sources64/x13as-$(X13AS_VERSION)

clean-x13as-bin:
	rm -rf lib64/x13as

clean-x13as-all: clean-x13as-tar clean-x13as-src clean-x13as-bin
