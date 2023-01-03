#!/usr/bin/env bash

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

set -ex

ROOTDIR=$(pwd)/..

# Set the GCC version
GCC_VERSION=12

# Set the compilers
CC=gcc-$GCC_VERSION
CXX=g++-$GCC_VERSION

# Set the number of threads
NTHREADS=$(nproc)

# Set dependency directory
LIB64="$ROOTDIR"/macOS/deps/lib64


##
## Find Dynare Version
##
DATE=$(date +%Y-%m-%d-%H%M)
DATELONG=$(date '+%d %B %Y')
if [[ -d ../.git/ ]]; then
    SHA=$(git rev-parse HEAD)
    SHASHORT=$(git rev-parse --short HEAD)
fi

if [[ -z $VERSION ]]; then
    VERSION=$(grep '^AC_INIT(' ../configure.ac | sed 's/AC_INIT(\[dynare\], \[\(.*\)\])/\1/')
    if [[ -d ../.git/ ]]; then
        VERSION="$VERSION"-"$SHASHORT"
    fi
fi

# Install location must not be too long for gcc.
# Otherwise, the headers of the compiled libraries cannot be modified
# obliging recompilation on the user's system.
# If VERSION is not a official release number, then do some magic
# to get something not too long, still more or less unique.
if [[ "$VERSION" =~ ^[0-9\.]+$ ]]; then
    LOCATION=$VERSION
else
    # Get the first component, truncate it to 5 characters, and add the date
    LOCATION=$(echo "$VERSION" | cut -f1 -d"-" | cut -c 1-5)-"$DATE"
fi

## Hack for statically linking libquadmath, similar to the one used in
## deps/Makefile for several libraries (there is no -static-libquadmath flag,
## see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=46539).
##
## NB: The hack done for Windows (see m4/ax_mexopts.m4) does not work for two reasons:
## - the macOS linker is different from GNU ld and does not have the equivalent of -Bstatic/-Bdynamic
## - libgfortran.spec does not include --as-needed on macOS, hence it will link the library anyways
## Also, it does not seem possible to override libgfortran.spec with the --specs option.
QUADMATH_DIR=$(mktemp -d)
ln -s /usr/local/opt/gcc/lib/gcc/$GCC_VERSION/libquadmath.a $QUADMATH_DIR

##
## Compile Dynare doc, preprocessor, mex for MATLAB < 2018a
##
## NB: In Homebrew, -static-libgfortran is implied by -static-libgcc (see “gfortran -dumpspecs”)
## NB2: We use the hack for libquadmath in LDFLAGS
cd "$ROOTDIR"
[[ -f configure ]] || autoreconf -si
./configure \
  PACKAGE_VERSION="$VERSION" \
  PACKAGE_STRING="dynare $VERSION" \
  CC=$CC \
  CXX=$CXX \
  CPPFLAGS=-I/usr/local/include \
  LDFLAGS="-static-libgcc -L$QUADMATH_DIR" \
  LEX=/usr/local/opt/flex/bin/flex \
  YACC=/usr/local/opt/bison/bin/bison \
  --with-gsl="$LIB64"/gsl \
  --with-matio="$LIB64"/matio \
  --with-slicot="$LIB64"/Slicot/with-underscore \
  --disable-octave \
  --with-matlab=/Applications/MATLAB_R2016b.app
if [[ -z $CI ]]; then
    # If not in Gitlab CI, clean the source and build the doc
    make clean
    make -j"$NTHREADS" pdf html
fi
make -j"$NTHREADS"


##
## Create package
##
NAME=dynare-"$VERSION"
PKGFILES="$ROOTDIR"/macOS/pkg/"$NAME"
mkdir -p \
      "$PKGFILES"/preprocessor \
      "$PKGFILES"/mex/matlab/maci64-8.3-9.3 \
      "$PKGFILES"/mex/matlab/maci64-9.4-9.13 \
      "$PKGFILES"/doc \
      "$PKGFILES"/scripts \
      "$PKGFILES"/contrib/ms-sbvar/TZcode

if [[ $VERSION == *-unstable* ]]; then
    echo "$SHA"                                                    > "$PKGFILES"/sha.txt
fi
cp -p  "$ROOTDIR"/NEWS.md                                            "$PKGFILES"
cp -p  "$ROOTDIR"/COPYING                                            "$PKGFILES"
cp -p  "$ROOTDIR"/VERSION                                            "$PKGFILES"
cp -p  "$ROOTDIR"/license.txt                                        "$PKGFILES"

cp -pr "$ROOTDIR"/matlab                                             "$PKGFILES"
cp -pr "$ROOTDIR"/examples                                           "$PKGFILES"

cp -p  "$ROOTDIR"/preprocessor/src/dynare-preprocessor               "$PKGFILES"/preprocessor

# Recreate backward-compatibility symlink
rm -f "$ROOTDIR"/matlab/preprocessor64/dynare_m
ln -sf ../../preprocessor/dynare-preprocessor                        "$PKGFILES"/matlab/preprocessor64/dynare_m

cp -L  "$ROOTDIR"/mex/matlab/*                                       "$PKGFILES"/mex/matlab/maci64-8.3-9.3

cp -p  "$ROOTDIR"/scripts/dynare.el                                  "$PKGFILES"/scripts
cp -pr "$ROOTDIR"/contrib/ms-sbvar/TZcode/MatlabFiles                "$PKGFILES"/contrib/ms-sbvar/TZcode
cp -pr "$ROOTDIR"/contrib/jsonlab                                    "$PKGFILES"/contrib

cp     "$ROOTDIR"/doc/*.pdf                                          "$PKGFILES"/doc
cp     "$ROOTDIR"/doc/gsa/gsa.pdf                                    "$PKGFILES"/doc
cp     "$ROOTDIR"/doc/parallel/parallel.pdf                          "$PKGFILES"/doc
cp     "$ROOTDIR"/doc/dseries-and-reporting/dseriesReporting.pdf     "$PKGFILES"/doc
cp     "$ROOTDIR"/preprocessor/doc/preprocessor/preprocessor.pdf     "$PKGFILES"/doc
cp     "$ROOTDIR"/preprocessor/doc/macroprocessor/macroprocessor.pdf "$PKGFILES"/doc
cp     "$ROOTDIR"/doc/manual/build/latex/dynare-manual.pdf           "$PKGFILES"/doc
cp -r  "$ROOTDIR"/doc/manual/build/html                              "$PKGFILES"/doc/dynare-manual.html

mkdir -p                                                             "$PKGFILES"/matlab/modules/dseries/externals/x13/macOS/64
cp -p  "$ROOTDIR"/macOS/deps/lib64/x13as/x13as                       "$PKGFILES"/matlab/modules/dseries/externals/x13/macOS/64


##
## Create mex for MATLAB ≥ 2018a
##
cd "$ROOTDIR"/mex/build/matlab
make clean
./configure \
  PACKAGE_VERSION="$VERSION" \
  PACKAGE_STRING="dynare $VERSION" \
  CC=$CC \
  CXX=$CXX \
  CPPFLAGS=-I/usr/local/include \
  LDFLAGS="-static-libgcc -L$QUADMATH_DIR" \
  --with-gsl="$LIB64"/gsl \
  --with-matio="$LIB64"/matio \
  --with-slicot="$LIB64"/Slicot/with-underscore \
  --with-matlab=/Applications/MATLAB_R2019b.app
make -j"$NTHREADS"
cp -L  "$ROOTDIR"/mex/matlab/*                                       "$PKGFILES"/mex/matlab/maci64-9.4-9.13


##
## Make package
##
cd "$ROOTDIR"/macOS/pkg

# Dynare option
pkgbuild --root "$PKGFILES" --identifier org.dynare --version "$VERSION" --install-location /Applications/Dynare/"$LOCATION" "$NAME".pkg

# GCC option
# Create dummy payload for GCC package; otherwise the size is displayed as 0 bytes in the installer
dd if=/dev/zero of="$ROOTDIR"/macOS/brewfiles/dummy  bs=1m  count=800
pkgbuild --root "$ROOTDIR"/macOS/brewfiles --identifier org.dynare.gcc --version "$VERSION" --scripts "$ROOTDIR"/macOS/scripts --install-location /Applications/Dynare/"$LOCATION" "$NAME"-gcc.pkg

# Replace variables in displayed files
sed "s/VERSION_READ/$VERSION/g" "$ROOTDIR"/macOS/distribution_template.xml > distribution_tmp.xml
sed "s/VERSION_NO_SPACE/$VERSION/g" distribution_tmp.xml > distribution.xml
sed "s/GCC_BINARY/$CC/g" "$ROOTDIR"/macOS/welcome_template.html > "$ROOTDIR"/macOS/welcome.html
sed "s/VERSION_NO_SPACE/$VERSION/g" "$ROOTDIR"/macOS/welcome.html > "$ROOTDIR"/macOS/welcome_tmp.html
sed "s/DATE/$DATELONG/g" "$ROOTDIR"/macOS/welcome_tmp.html > "$ROOTDIR"/macOS/welcome.html

# Create installer
productbuild --distribution distribution.xml --resources "$ROOTDIR"/macOS --package-path ./"$NAME".pkg "$NAME"-new.pkg

# cleanup
rm -f ./*.xml
rm -rf "$PKGFILES"
rm -f "$NAME"-gcc.pkg
rm -f "$ROOTDIR"/macOS/brewfiles/dummy
rm -f "$ROOTDIR"/macOS/welcome.html
rm -f "$ROOTDIR"/macOS/welcome_tmp.html

# Final pkg
mv "$NAME"-new.pkg "$NAME".pkg
