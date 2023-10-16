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
#exec > >(tee build-logfile.log) 2>&1 # uncomment for debugging

ROOTDIR=$(pwd)/..
##
## Set settings based on architecture
##
path_remove ()  { export $1="`echo -n ${!1} | awk -v RS=: -v ORS=: '$1 != "'$2'"' | sed 's/:$//'`"; }
path_prepend () { path_remove $1 $2; export $1="$2:${!1}"; }
PKG_ARCH=$1
if [[ $PKG_ARCH == arm64 ]]; then
    BREWDIR=/opt/homebrew
    # Make sure /opt/homebrew/bin is set first in PATH (as it might come last)
    path_prepend PATH /opt/homebrew/bin
    MATLAB_ARCH=maca64
    # arm64 MATLAB is only available starting with R2023b, no need to distinguish versions
    MATLAB_PATH_BASE=/Applications/"$PKG_ARCH"/MATLAB_R2023b.app
else
    BREWDIR=/usr/local
    # Remove /opt/homebrew/bin from PATH, so it does not intervene with the x86_64 compilations
    path_remove PATH /opt/homebrew/bin
    MATLAB_ARCH=maci64
    # On x86_64 we need to differentiate between older (BASE) and newer (NEW) MATLAB versions due to ABI breaks
    MATLAB_PATH_BASE=/Applications/MATLAB_R2016b.app
    MATLAB_PATH_NEW=/Applications/"$PKG_ARCH"/MATLAB_R2023b.app
fi

# Append texbin to PATH to access latexmk and friends
path_prepend PATH /Library/TeX/texbin

# Set the GCC version
GCC_VERSION=13

# Set the compilers
CC=gcc-$GCC_VERSION
CXX=g++-$GCC_VERSION

# Set the number of threads
NTHREADS=$(sysctl -n hw.ncpu)

# Set dependency directory
LIB64="$ROOTDIR"/macOS/deps/"$PKG_ARCH"/lib64


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
# Add architecture to LOCATION and VERSION
VERSION="$VERSION"-"$PKG_ARCH"
LOCATION="$LOCATION"-"$PKG_ARCH"

## Hack for statically linking libquadmath, similar to the one used in
## deps/Makefile for several libraries (there is no -static-libquadmath flag,
## see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=46539).
##
## NB: The hack done for Windows (see m4/ax_mexopts.m4) does not work for two reasons:
## - the macOS linker is different from GNU ld and does not have the equivalent of -Bstatic/-Bdynamic
## - libgfortran.spec does not include --as-needed on macOS, hence it will link the library anyways
## Also, it does not seem possible to override libgfortran.spec with the --specs option.
QUADMATH_DIR=$(mktemp -d)
ln -s $BREWDIR/opt/gcc/lib/gcc/$GCC_VERSION/libquadmath.a $QUADMATH_DIR

##
## Compile Dynare doc, dynare++, preprocessor, mex for MATLAB < 2018a (x86_64) or MATLAB = R2023b (arm64)
##
## NB: In Homebrew, -static-libgfortran is implied by -static-libgcc (see “gfortran -dumpspecs”)
## NB2: We use the hack for libquadmath in LDFLAGS
## NB3: The -Wl,-ld_classic flag is a workaround for a bug in XCode 15
cd "$ROOTDIR"
[[ -f configure ]] || arch -$PKG_ARCH autoreconf -si
arch -$PKG_ARCH ./configure \
  PACKAGE_VERSION="$VERSION" \
  PACKAGE_STRING="dynare $VERSION" \
  CC=$CC \
  CXX=$CXX \
  CPPFLAGS=-I$BREWDIR/include \
  LDFLAGS="-static-libgcc -L$QUADMATH_DIR -Wl,-ld_classic" \
  LEX=$BREWDIR/opt/flex/bin/flex \
  YACC=$BREWDIR/opt/bison/bin/bison \
  --with-gsl="$LIB64"/gsl \
  --with-matio="$LIB64"/matio \
  --with-slicot="$LIB64"/slicot/with-underscore \
  --disable-octave \
  --with-matlab="$MATLAB_PATH_BASE"
if [[ -z $CI ]]; then
    # If not in Gitlab CI, clean the source and build the doc
    arch -$PKG_ARCH make clean
    arch -$PKG_ARCH make -j"$NTHREADS" pdf html
fi
arch -$PKG_ARCH make -j"$NTHREADS"


##
## Create package
##
NAME=dynare-"$VERSION"
PKGFILES="$ROOTDIR"/macOS/pkg/"$NAME"
mkdir -p \
      "$PKGFILES"/preprocessor \
      "$PKGFILES"/doc/dynare++ \
      "$PKGFILES"/dynare++ \
      "$PKGFILES"/scripts \
      "$PKGFILES"/contrib/ms-sbvar/TZcode
if [[ $PKG_ARCH == x86_64 ]]; then
    mkdir -p "$PKGFILES"/mex/matlab/"$MATLAB_ARCH"-8.3-9.3 \
             "$PKGFILES"/mex/matlab/"$MATLAB_ARCH"-9.4-23.2
else
    mkdir -p "$PKGFILES"/mex/matlab/"$MATLAB_ARCH"-23.2
fi      

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

if [[ $PKG_ARCH == x86_64 ]]; then    
    cp -L  "$ROOTDIR"/mex/matlab/*                                   "$PKGFILES"/mex/matlab/"$MATLAB_ARCH"-8.3-9.3
else
    cp -L  "$ROOTDIR"/mex/matlab/*                                   "$PKGFILES"/mex/matlab/"$MATLAB_ARCH"-23.2
fi

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

cp     "$ROOTDIR"/dynare++/doc/*.pdf                                 "$PKGFILES"/doc/dynare++

cp     "$ROOTDIR"/dynare++/src/dynare++                              "$PKGFILES"/dynare++
cp     "$ROOTDIR"/dynare++/dynare_simul/dynare_simul.m               "$PKGFILES"/dynare++

mkdir -p                                                             "$PKGFILES"/matlab/modules/dseries/externals/x13/macOS/64
cp -p  "$ROOTDIR"/macOS/deps/$PKG_ARCH/lib64/x13as/x13as             "$PKGFILES"/matlab/modules/dseries/externals/x13/macOS/64


##
## Create mex for MATLAB ≥ 2018a (only for x86_64)
##
if [[ $PKG_ARCH == x86_64 ]]; then
    cd "$ROOTDIR"/mex/build/matlab
    arch -$PKG_ARCH make clean
    arch -$PKG_ARCH ./configure \
      PACKAGE_VERSION="$VERSION" \
      PACKAGE_STRING="dynare $VERSION" \
      CC=$CC \
      CXX=$CXX \
      CPPFLAGS=-I$BREWDIR/include \
      LDFLAGS="-static-libgcc -L$QUADMATH_DIR -Wl,-ld_classic" \
      --with-gsl="$LIB64"/gsl \
      --with-matio="$LIB64"/matio \
      --with-slicot="$LIB64"/slicot/with-underscore \
      --with-matlab="$MATLAB_PATH_NEW"
    arch -$PKG_ARCH make -j"$NTHREADS"
    cp -L  "$ROOTDIR"/mex/matlab/*                                   "$PKGFILES"/mex/matlab/$MATLAB_ARCH-9.4-23.2
fi

##
## Make package
##
cd "$ROOTDIR"/macOS/pkg

# Dynare option
arch -$PKG_ARCH pkgbuild --root "$PKGFILES" --identifier org.dynare."$VERSION" --version "$VERSION" --install-location /Applications/Dynare/"$LOCATION" "$NAME".pkg

# Create distribution.xml by replacing variables in distribution_template.xml
sed -e "s/VERSION_NO_SPACE/$VERSION/g" \
    -e "s/LOCATION/$LOCATION/g" \
    "$ROOTDIR"/macOS/distribution_template.xml > distribution.xml

# Create welcome.html by replacing variables in welcome_template.html
sed -e "s/VERSION_NO_SPACE/$VERSION/g" \
    -e "s/DATE/$DATELONG/g" \
    -e "s/GCC_VERSION/$GCC_VERSION/g" \
    "$ROOTDIR"/macOS/welcome_template.html > "$ROOTDIR"/macOS/welcome.html

# Create conclusion.html by replacing variables in conclusion_template.html
sed -e "s/GCC_VERSION/$GCC_VERSION/g" \
    "$ROOTDIR"/macOS/conclusion_template.html > "$ROOTDIR"/macOS/conclusion.html

# Create installer
arch -$PKG_ARCH productbuild --distribution distribution.xml --resources "$ROOTDIR"/macOS --package-path ./"$NAME".pkg "$NAME"-productbuild.pkg

# Cleanup
rm -f ./distribution.xml
rm -rf "$PKGFILES"
rm -f "$ROOTDIR"/macOS/welcome.html
rm -f "$ROOTDIR"/macOS/conclusion.html

# Final pkg
mv "$NAME"-productbuild.pkg "$NAME".pkg
