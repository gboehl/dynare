# Copyright © 2017-2019 Dynare Team
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

include versions.mk

ROOT_PATH = $(realpath .)

WGET_OPTIONS := --no-verbose --no-use-server-timestamps --retry-connrefused --retry-on-host-error

.PHONY: clean-openblas-src clean-openblas-tar clean-openblas-all clean-libopenblas \
	clean-boost-src clean-boost-tar clean-boost-all clean-libboost \
	clean-gsl-src clean-gsl-tar clean-gsl-all clean-libgsl\
	clean-slicot-src clean-slicot-tar clean-slicot-all clean-libslicot \
	clean-matio-src clean-matio-tar clean-matio-all clean-libmatio\
	clean-zlib-src clean-zlib-tar clean-zlib-src clean-libzlib\
	clean-src clean-all clean-lib clean-tar \
	build build-openblas build-slicot build-matio build-boost build-gsl \
	octave matlab mingw

all: build octave matlab mingw

build: build-openblas build-slicot build-matio build-boost build-gsl

clean-lib: clean-libopenblas clean-libgsl clean-libzlib clean-libmatio clean-libslicot clean-libboost clean-matlab clean-octave clean-mingw

clean-src: clean-openblas-src clean-boost-src clean-gsl-src clean-matio-src clean-slicot-src clean-zlib-src

clean-tar: clean-openblas-tar clean-boost-tar clean-gsl-tar clean-matio-tar clean-slicot-tar clean-zlib-tar

clean-all: clean-lib clean-src clean-tar

#
# OpenBLAS library
#

tarballs/openblas-${OPENBLAS_VERSION}.tar.gz:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://github.com/xianyi/OpenBLAS/archive/v${OPENBLAS_VERSION}.tar.gz

sources32/OpenBLAS-${OPENBLAS_VERSION}: tarballs/openblas-${OPENBLAS_VERSION}.tar.gz
	mkdir -p sources32
	tar xf $< --directory sources32
	patch $@/Makefile.rule < patches/openblas-w32.patch

sources64/OpenBLAS-${OPENBLAS_VERSION}: tarballs/openblas-${OPENBLAS_VERSION}.tar.gz
	mkdir -p sources64
	tar xf $< --directory sources64
	patch $@/Makefile.rule < patches/openblas-w64.patch

lib32/OpenBLAS/libopenblas.a: sources32/OpenBLAS-${OPENBLAS_VERSION}
	make -C $<
	i686-w64-mingw32-strip --strip-debug $</libopenblasp-r${OPENBLAS_VERSION}.a
	mkdir -p $(dir $@)
	cp $</libopenblasp-r${OPENBLAS_VERSION}.a $@

lib64/OpenBLAS/libopenblas.a: sources64/OpenBLAS-${OPENBLAS_VERSION}
	make -C $<
	x86_64-w64-mingw32-strip --strip-debug $</libopenblasp-r${OPENBLAS_VERSION}.a
	mkdir -p $(dir $@)
	cp $</libopenblasp-r${OPENBLAS_VERSION}.a $@

build-openblas: lib32/OpenBLAS/libopenblas.a lib64/OpenBLAS/libopenblas.a

clean-openblas-tar:
	rm -f tarballs/openblas-${OPENBLAS_VERSION}.tar.gz

clean-openblas-src:
	rm -rf sources32/OpenBLAS-${OPENBLAS_VERSION}
	rm -rf sources64/OpenBLAS-${OPENBLAS_VERSION}

clean-libopenblas:
	rm -rf lib32/OpenBLAS
	rm -rf lib64/OpenBLAS

clean-openblas-all: clean-openblas-src clean-openblas-tar clean-libopenblas

#
# Boost library
#

tarballs/boost_${BOOST_VERSION}.tar.bz2:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://dl.bintray.com/boostorg/release/$(subst _,.,${BOOST_VERSION})/source/$(notdir $@)

sources32/boost_${BOOST_VERSION}: tarballs/boost_${BOOST_VERSION}.tar.bz2
	mkdir -p sources32
	tar xf $< --directory sources32
	echo "using gcc : mingw : i686-w64-mingw32-g++ ;" > $@/user-config.jam

sources64/boost_${BOOST_VERSION}: tarballs/boost_${BOOST_VERSION}.tar.bz2
	mkdir -p sources64
	tar xf $< --directory sources64
	echo "using gcc : mingw : x86_64-w64-mingw32-g++ ;" > $@/user-config.jam

lib32/Boost: sources32/boost_${BOOST_VERSION}
	rm -rf lib32/Boost
	cd $< && ./bootstrap.sh --with-libraries=system,filesystem --prefix=${ROOT_PATH}/lib32/Boost && ./b2 -a -q --user-config=user-config.jam toolset=gcc-mingw target-os=windows address-model=32 variant=release install

lib64/Boost: sources64/boost_${BOOST_VERSION}
	rm -rf lib64/Boost
	cd $< && ./bootstrap.sh --with-libraries=system,filesystem --prefix=${ROOT_PATH}/lib64/Boost && ./b2 -a -q --user-config=user-config.jam toolset=gcc-mingw target-os=windows address-model=64 variant=release install

build-boost: lib32/Boost lib64/Boost

clean-boost-tar:
	rm -f tarballs/boost_${BOOST_VERSION}.tar.bz2

clean-boost-src:
	rm -rf sources32/boost_${BOOST_VERSION}
	rm -rf sources64/boost_${BOOST_VERSION}

clean-libboost:
	rm -rf lib32/Boost
	rm -rf lib64/Boost

clean-boost-all: clean-boost-src clean-boost-tar clean-libboost

#
# Gsl
#

tarballs/gsl-${GSL_VERSION}.tar.gz:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://ftp.gnu.org/gnu/gsl/$(notdir $@)

sources32/gsl-${GSL_VERSION}: tarballs/gsl-${GSL_VERSION}.tar.gz
	mkdir -p sources32
	tar xf $< --directory sources32

sources64/gsl-${GSL_VERSION}: tarballs/gsl-${GSL_VERSION}.tar.gz
	mkdir -p sources64
	tar xf $< --directory sources64

lib32/Gsl/lib/libgsl.a: sources32/gsl-${GSL_VERSION}
	rm -rf lib32/Gsl
	cd $< && ./configure --host=i686-w64-mingw32 --prefix=${ROOT_PATH}/lib32/Gsl --disable-shared --enable-static && make && make install

lib64/Gsl/lib/libgsl.a: sources64/gsl-${GSL_VERSION}
	rm -rf lib64/Gsl
	cd $< && ./configure --host=x86_64-w64-mingw32 --prefix=${ROOT_PATH}/lib64/Gsl --disable-shared --enable-static && make && make install

build-gsl: lib32/Gsl/lib/libgsl.a lib64/Gsl/lib/libgsl.a

clean-gsl-tar:
	rm -f tarballs/gsl-${GSL_VERSION}.tar.gz

clean-gsl-src:
	rm -rf sources32/gsl-${GSL_VERSION}
	rm -rf sources64/gsl-${GSL_VERSION}

clean-libgsl:
	rm -rf lib32/Gsl
	rm -rf lib64/Gsl

clean-gsl-all: clean-gsl-src clean-gsl-tar clean-libgsl

#
# matIO
#

tarballs/matio-${MATIO_VERSION}.tar.gz:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://sourceforge.net/projects/matio/files/matio/${MATIO_VERSION}/$(notdir $@)/download

sources32/matio-${MATIO_VERSION}: tarballs/matio-${MATIO_VERSION}.tar.gz
	mkdir -p sources32
	tar xf $< --directory sources32
	touch $@

sources64/matio-${MATIO_VERSION}: tarballs/matio-${MATIO_VERSION}.tar.gz
	mkdir -p sources64
	tar xf $< --directory sources64
	touch $@

# The "ac_cv_va_copy=C99" argument is a workaround for https://github.com/tbeu/matio/issues/78

lib32/matIO/lib/libmatio.a: sources32/matio-${MATIO_VERSION} lib32/Zlib/lib/libz.a
	rm -rf lib32/matIO
	cd $< && ./configure --host=i686-w64-mingw32 --disable-shared --with-zlib=${ROOT_PATH}/lib32/Zlib --prefix=${ROOT_PATH}/lib32/matIO ac_cv_va_copy=C99 && make install
	ln -sf ${ROOT_PATH}/lib32/Zlib/include/zconf.h ${ROOT_PATH}/lib32/matIO/include/zconf.h
	ln -sf ${ROOT_PATH}/lib32/Zlib/include/zlib.h ${ROOT_PATH}/lib32/matIO/include/zlib.h
	ln -sf ${ROOT_PATH}/lib32/Zlib/lib/libz.a ${ROOT_PATH}/lib32/matIO/lib/libz.a

lib64/matIO/lib/libmatio.a: sources64/matio-${MATIO_VERSION} lib64/Zlib/lib/libz.a
	rm -rf lib64/matIO
	cd $< && ./configure --host=x86_64-w64-mingw32 --disable-shared --with-zlib=${ROOT_PATH}/lib64/Zlib --prefix=${ROOT_PATH}/lib64/matIO ac_cv_va_copy=C99 && make install
	ln -sf ${ROOT_PATH}/lib64/Zlib/include/zconf.h ${ROOT_PATH}/lib64/matIO/include/zconf.h
	ln -sf ${ROOT_PATH}/lib64/Zlib/include/zlib.h ${ROOT_PATH}/lib64/matIO/include/zlib.h
	ln -sf ${ROOT_PATH}/lib64/Zlib/lib/libz.a ${ROOT_PATH}/lib64/matIO/lib/libz.a

build-matio: build-zlib lib32/matIO/lib/libmatio.a lib64/matIO/lib/libmatio.a

clean-matio-tar:
	rm -f tarballs/matio-${MATIO_VERSION}.tar.gz

clean-matio-src:
	rm -rf sources32/matio-${MATIO_VERSION}
	rm -rf sources64/matio-${MATIO_VERSION}

clean-libmatio:
	rm -rf lib32/matIO
	rm -rf lib64/matIO

clean-matio-all: clean-matio-src clean-matio-tar clean-libmatio

#
# Slicot
#

tarballs/slicot-$(SLICOT_VERSION).tar.gz:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://deb.debian.org/debian/pool/main/s/slicot/slicot_$(SLICOT_VERSION).orig.tar.gz

sources32/slicot-$(SLICOT_VERSION)-without-underscore: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1
	patch $@/make.inc < patches/slicot-32-without-underscore.patch

sources32/slicot-$(SLICOT_VERSION)-with-underscore: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1
	patch $@/make.inc < patches/slicot-32-with-underscore.patch

sources64/slicot-$(SLICOT_VERSION)-with-32bit-integer: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1
	patch $@/make.inc < patches/slicot-64-with-32bit-integer.patch

sources64/slicot-$(SLICOT_VERSION)-with-64bit-integer: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1
	patch $@/make.inc < patches/slicot-64-with-64bit-integer.patch

sources64/slicot-$(SLICOT_VERSION)-with-32bit-integer-and-underscore: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1
	patch $@/make.inc < patches/slicot-64-with-32bit-integer-and-underscore.patch

sources64/slicot-$(SLICOT_VERSION)-with-64bit-integer-and-underscore: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1
	patch $@/make.inc < patches/slicot-64-with-64bit-integer-and-underscore.patch

lib32/Slicot/without-underscore/lib/libslicot_pic.a: sources32/slicot-$(SLICOT_VERSION)-without-underscore
	make -C $< lib
	i686-w64-mingw32-strip --strip-debug $</libslicot_pic.a
	mkdir -p $(dir $@)
	cp $</libslicot_pic.a $@

lib32/Slicot/with-underscore/lib/libslicot_pic.a: sources32/slicot-$(SLICOT_VERSION)-with-underscore
	make -C $< lib
	i686-w64-mingw32-strip --strip-debug $</libslicot_pic.a
	mkdir -p $(dir $@)
	cp $</libslicot_pic.a $@

lib64/Slicot/without-underscore/lib/libslicot_pic.a: sources64/slicot-$(SLICOT_VERSION)-with-32bit-integer
	make -C $< lib
	x86_64-w64-mingw32-strip --strip-debug $</libslicot_pic.a
	mkdir -p $(dir $@)
	cp $</libslicot_pic.a $@

lib64/Slicot/without-underscore/lib/libslicot64_pic.a: sources64/slicot-$(SLICOT_VERSION)-with-64bit-integer
	make -C $< lib
	x86_64-w64-mingw32-strip --strip-debug $</libslicot64_pic.a
	mkdir -p $(dir $@)
	cp $</libslicot64_pic.a $@

lib64/Slicot/with-underscore/lib/libslicot_pic.a: sources64/slicot-$(SLICOT_VERSION)-with-32bit-integer-and-underscore
	make -C $< lib
	x86_64-w64-mingw32-strip --strip-debug $</libslicot_pic.a
	mkdir -p $(dir $@)
	cp $</libslicot_pic.a $@

lib64/Slicot/with-underscore/lib/libslicot64_pic.a: sources64/slicot-$(SLICOT_VERSION)-with-64bit-integer-and-underscore
	make -C $< lib
	x86_64-w64-mingw32-strip --strip-debug $</libslicot64_pic.a
	mkdir -p $(dir $@)
	cp $</libslicot64_pic.a $@

build-slicot: lib32/Slicot/without-underscore/lib/libslicot_pic.a \
	lib32/Slicot/with-underscore/lib/libslicot_pic.a \
	lib64/Slicot/without-underscore/lib/libslicot_pic.a \
	lib64/Slicot/without-underscore/lib/libslicot64_pic.a \
	lib64/Slicot/with-underscore/lib/libslicot_pic.a \
	lib64/Slicot/with-underscore/lib/libslicot64_pic.a

clean-slicot-tar:
	rm -f tarballs/slicot-$(SLICOT_VERSION).tar.gz

clean-slicot-src:
	rm -rf sources32/slicot-$(SLICOT_VERSION)-with-underscore
	rm -rf sources32/slicot-$(SLICOT_VERSION)-without-underscore
	rm -rf sources64/slicot-$(SLICOT_VERSION)-with-32bit-integer
	rm -rf sources64/slicot-$(SLICOT_VERSION)-with-64bit-integer
	rm -rf sources64/slicot-$(SLICOT_VERSION)-with-32bit-integer-and-underscore
	rm -rf sources64/slicot-$(SLICOT_VERSION)-with-64bit-integer-and-underscore

clean-libslicot:
	rm -rf lib64/Slicot
	rm -rf lib32/Slicot

clean-slicot-all: clean-slicot-src clean-slicot-tar clean-libslicot

#
# Zlib
#

tarballs/zlib-${ZLIB_VERSION}.tar.xz:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://zlib.net/$(notdir $@)

sources32/zlib-${ZLIB_VERSION}: tarballs/zlib-${ZLIB_VERSION}.tar.xz
	mkdir -p sources32
	tar xf $< --directory sources32
	touch $@

sources64/zlib-${ZLIB_VERSION}: tarballs/zlib-${ZLIB_VERSION}.tar.xz
	mkdir -p sources64
	tar xf $< --directory sources64
	touch $@

lib32/Zlib/lib/libz.a: sources32/zlib-${ZLIB_VERSION}
	rm -rf lib32/Zlib
	cd $< && CC=i686-w64-mingw32-gcc ./configure --static --prefix=${ROOT_PATH}/lib32/Zlib && make install

lib64/Zlib/lib/libz.a: sources64/zlib-${ZLIB_VERSION}
	rm -rf lib64/Zlib
	cd $< && CC=x86_64-w64-mingw32-gcc ./configure --static --prefix=${ROOT_PATH}/lib64/Zlib && make install

build-zlib: lib32/Zlib/lib/libz.a lib64/Zlib/lib/libz.a

clean-zlib-tar:
	rm -f tarballs/zlib-${ZLIB_VERSION}.tar.xz

clean-zlib-src:
	rm -rf sources32/zlib-${ZLIB_VERSION}
	rm -rf sources64/zlib-${ZLIB_VERSION}

clean-libzlib:
	rm -rf lib32/Zlib
	rm -rf lib64/Zlib

clean-zlib-all: clean-zlib-src clean-zlib-tar clean-libzlib

#
# Octave
#

tarballs/octave-${OCTAVE_VERSION}-w32${OCTAVE_W32_BUILD}.7z:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://ftp.gnu.org/gnu/octave/windows/$(notdir $@)

tarballs/octave-${OCTAVE_VERSION}-w64${OCTAVE_W64_BUILD}.7z:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://ftp.gnu.org/gnu/octave/windows/$(notdir $@)

lib32/octave/mingw32/octave-cli.exe: tarballs/octave-${OCTAVE_VERSION}-w32${OCTAVE_W32_BUILD}.7z
	rm -rf lib32/octave
	7zr x -olib32 $< > /dev/null
	mv lib32/octave-${OCTAVE_VERSION}-w32 lib32/octave
	# Clean up MinGW stuff that interferes with our cross-compiler
	cd lib32/octave/mingw32/lib && rm -f libadvapi32.a libkernel32.a libmingw32.a libmingwex.a libmoldname.a libmsvcrt.a libpthread.dll.a libshell32.a libuser32.a
	cd lib32/octave/mingw32/include && rm -f $$(grep -l "This file is part of the mingw-w64 runtime package" *.h)
	touch $@

lib64/octave/mingw64/octave-cli.exe: tarballs/octave-${OCTAVE_VERSION}-w64${OCTAVE_W64_BUILD}.7z
	rm -rf lib64/octave
	7zr x -olib64 $< > /dev/null
	mv lib64/octave-${OCTAVE_VERSION}-w64 lib64/octave
	# Clean up MinGW stuff that interferes with our cross-compiler
	cd lib64/octave/mingw64/lib && rm -f libadvapi32.a libkernel32.a libmingw32.a libmingwex.a libmoldname.a libmsvcrt.a libpthread.dll.a libshell32.a libuser32.a
	cd lib64/octave/mingw64/include && rm -f $$(grep -l "This file is part of the mingw-w64 runtime package" *.h)
	touch $@

octave: lib32/octave/mingw32/octave-cli.exe lib64/octave/mingw64/octave-cli.exe

clean-octave:
	rm -rf lib32/octave
	rm -rf lib64/octave
	rm -f tarballs/octave-${OCTAVE_VERSION}-w64${OCTAVE_W64_BUILD}.7z
	rm -f tarballs/octave-${OCTAVE_VERSION}-w32${OCTAVE_W32_BUILD}.7z

#
# MATLAB
#

tarballs/matlab%.tar.xz:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@.gpg https://www.dynare.org/windows-pkg-build/$(notdir $@).gpg
	gpg --output $@ --decrypt $@.gpg
	rm $@.gpg

lib32/matlab/R2009b/bin/win32/MATLAB.exe: tarballs/matlab32-$(MATLAB32_VERSION).tar.xz
	rm -rf lib32/matlab
	mkdir -p lib32/matlab
	tar xf $< --directory lib32/matlab --strip-components=1
	touch $@

lib64/matlab/R2009b/bin/win64/MATLAB.exe: tarballs/matlab64-$(MATLAB64_VERSION).tar.xz
	rm -rf lib64/matlab
	mkdir -p lib64/matlab
	tar xf $< --directory lib64/matlab --strip-components=1
	touch $@

matlab: lib32/matlab/R2009b/bin/win32/MATLAB.exe lib64/matlab/R2009b/bin/win64/MATLAB.exe

clean-matlab:
	rm -f tarballs/matlab32.tar.xz
	rm -f tarballs/matlab64.tar.xz
	rm -rf lib32/matlab
	rm -rf lib64/matlab

#
# Embedded MinGW (for use_dll)
#

mingw: mingw32 mingw64

mingw32: tarballs/mingw-w64-i686-gcc-$(MINGW32_GCC_VERSION)-any.pkg.tar.xz tarballs/mingw-w64-i686-gmp-$(MINGW32_GMP_VERSION)-any.pkg.tar.xz tarballs/mingw-w64-i686-binutils-$(MINGW32_BINUTILS_VERSION)-any.pkg.tar.xz tarballs/mingw-w64-i686-headers-git-$(MINGW32_HEADERS_VERSION)-any.pkg.tar.xz tarballs/mingw-w64-i686-crt-git-$(MINGW32_CRT_VERSION)-any.pkg.tar.xz tarballs/mingw-w64-i686-winpthreads-git-$(MINGW32_WINPTHREADS_VERSION)-any.pkg.tar.xz tarballs/mingw-w64-i686-libwinpthread-git-$(MINGW32_LIBWINPTHREAD_VERSION)-any.pkg.tar.xz
	rm -rf $@
	for f in $^; do tar xf $$f; done
	touch $@

mingw64: tarballs/mingw-w64-x86_64-gcc-$(MINGW64_GCC_VERSION)-any.pkg.tar.xz tarballs/mingw-w64-x86_64-gmp-$(MINGW64_GMP_VERSION)-any.pkg.tar.xz tarballs/mingw-w64-x86_64-binutils-$(MINGW64_BINUTILS_VERSION)-any.pkg.tar.xz tarballs/mingw-w64-x86_64-headers-git-$(MINGW64_HEADERS_VERSION)-any.pkg.tar.xz tarballs/mingw-w64-x86_64-crt-git-$(MINGW64_CRT_VERSION)-any.pkg.tar.xz tarballs/mingw-w64-x86_64-winpthreads-git-$(MINGW64_WINPTHREADS_VERSION)-any.pkg.tar.xz tarballs/mingw-w64-x86_64-libwinpthread-git-$(MINGW64_LIBWINPTHREAD_VERSION)-any.pkg.tar.xz
	rm -rf $@
	for f in $^; do tar xf $$f; done
	touch $@

tarballs/mingw-w64-i686-%-any.pkg.tar.xz:
	wget $(WGET_OPTIONS) -O $@ http://repo.msys2.org/mingw/i686/$(notdir $@)

tarballs/mingw-w64-x86_64-%-any.pkg.tar.xz:
	wget $(WGET_OPTIONS) -O $@ http://repo.msys2.org/mingw/x86_64/$(notdir $@)

clean-mingw:
	rm -rf mingw32
	rm -rf mingw64
	rm -f tarballs/mingw-w64-i686-*.tar.xz
	rm -f tarballs/mingw-w64-x86_64-*.tar.xz
