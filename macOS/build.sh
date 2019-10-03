#!/usr/bin/env bash

# Copyright © 2019 Dynare Team
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

set -ex

# Set the compilers
CC=gcc-9
CXX=g++-9

# Set the number of threads
NTHREADS=$(nproc)

##
## Find Dynare Version
##
ROOTDIR=$(pwd)/..
if [[ -z $VERSION ]]; then
    VERSION=$(grep '^AC_INIT(' ../configure.ac | sed 's/AC_INIT(\[dynare\], \[\(.*\)\])/\1/')
    if [[ -d ../.git/ ]]; then
        SHA=$(git rev-parse --short HEAD)
        VERSION_READ="$VERSION-$SHA"
        VERSION=$VERSION-$(date +%Y-%m-%d-%H%M)-"$SHA"
    fi
fi

# Set dependency directory
LIB64="$ROOTDIR"/macOS/deps/lib64


##
## Compile Dynare
##
cd "$ROOTDIR"
[[ -f configure ]] || autoreconf -si
CC=$CC CXX=$CXX ./configure --with-matlab=/Applications/MATLAB_R2016b.app MATLAB_VERSION=R2016b --with-matio=/usr/local --with-gsl=/usr/local --with-slicot="$LIB64"/Slicot/with-underscore --disable-octave PACKAGE_VERSION="$VERSION" PACKAGE_STRING="dynare $VERSION"
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
mkdir -p "$PKGFILES"/mex/matlab/maci64-8.7-9.3
mkdir    "$PKGFILES"/mex/matlab/maci64-9.4-9.7
mkdir    "$PKGFILES"/mex/octave
mkdir -p "$PKGFILES"/doc/dynare++
mkdir    "$PKGFILES"/dynare++
mkdir    "$PKGFILES"/scripts
mkdir    "$PKGFILES"/contrib

cp -p  "$ROOTDIR"/NEWS                                               "$PKGFILES"
cp -p  "$ROOTDIR"/COPYING                                            "$PKGFILES"
cp -p  "$ROOTDIR"/VERSION                                            "$PKGFILES"
cp -p  "$ROOTDIR"/license.txt                                        "$PKGFILES"

cp -pr "$ROOTDIR"/matlab                                             "$PKGFILES"
cp -pr "$ROOTDIR"/examples                                           "$PKGFILES"

cp -L  "$ROOTDIR"/mex/matlab/*                                       "$PKGFILES"/mex/matlab/maci64-8.7-9.3

cp -p  "$ROOTDIR"/scripts/dynare.el                                  "$PKGFILES"/scripts
cp -pr "$ROOTDIR"/contrib/ms-sbvar                                   "$PKGFILES"/contrib
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

mkdir -p                                                             "$PKGFILES"/matlab/modules/dseries/externals/x13/osx/64
cp -p  "$ROOTDIR"/macOS/deps/lib64/x13as/x13as                       "$PKGFILES"/matlab/modules/dseries/externals/x13/osx/64


##
## Create mex for MATLAB le 2018a
##
cd "$ROOTDIR"/mex/build/matlab
make clean
CC=$CC CXX=$CXX ./configure --with-matlab=/Applications/MATLAB_R2019b.app MATLAB_VERSION=R2019b --with-matio=/usr/local --with-gsl=/usr/local --with-slicot="$LIB64"/Slicot/with-underscore PACKAGE_VERSION="$VERSION" PACKAGE_STRING="dynare $VERSION"
make -j"$NTHREADS"
cp -L  "$ROOTDIR"/mex/matlab/*                                       "$PKGFILES"/mex/matlab/maci64-9.4-9.7


##
## Create mex for Octave
##
cd "$ROOTDIR"/mex/build/octave
CC=$CC CXX=$CXX ./configure --with-matio=/usr/local --with-gsl=/usr/local --with-slicot="$LIB64"/Slicot/with-underscore LDFLAGS=-L/usr/local/lib PACKAGE_VERSION="$VERSION" PACKAGE_STRING="dynare $VERSION"
make -j"$NTHREADS"
cp -L  "$ROOTDIR"/mex/octave/*                                       "$PKGFILES"/mex/octave
echo -e "function v = supported_octave_version\nv=\"$(octave --eval "disp(OCTAVE_VERSION)")\";\nend" > "$PKGFILES"/matlab/supported_octave_version.m


##
## Make package
##
cd "$ROOTDIR"/macOS/pkg
pkgbuild --root "$PKGFILES" --identifier com.cepremap.dynare --version "$VERSION" --install-location /Applications/Dynare/"$VERSION" "$NAME".pkg
sed "s/VERSION_READ/$VERSION_READ/g" "$ROOTDIR"/macOS/distribution_template.xml > distribution_tmp.xml
sed "s/VERSION_NO_SPACE/$VERSION/g" distribution_tmp.xml > distribution.xml
ln -s "$ROOTDIR"/COPYING "$ROOTDIR"/macOS/
productbuild --distribution distribution.xml --resources "$ROOTDIR"/macOS --package-path ./"$NAME".pkg "$NAME"-new.pkg
rm -f ./*.xml
rm -rf "$PKGFILES"
rm "$ROOTDIR"/macOS/COPYING
mv "$NAME"-new.pkg "$NAME".pkg
