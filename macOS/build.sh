#!/usr/bin/env bash

set -ex

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


##
## Compile Dynare
##
cd "$ROOTDIR"
[[ -f configure ]] || autoreconf -si
CC=gcc-9 CXX=g++-9 ./configure --with-matlab=/Applications/MATLAB_R2016b.app MATLAB_VERSION=R2016b --with-matio=/usr/local --with-gsl=/usr/local --with-slicot=/usr/local --disable-octave PACKAGE_VERSION="$VERSION" PACKAGE_STRING="dynare $VERSION"
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


##
## Create mex for Matlab le 2018a
##
cd "$ROOTDIR"/mex/build/matlab
make clean
CC=gcc-9 CXX=g++-9 ./configure --with-matlab=/Applications/MATLAB_R2019b.app MATLAB_VERSION=R2019b --with-matio=/usr/local --with-gsl=/usr/local --with-slicot=/usr/local PACKAGE_VERSION="$VERSION" PACKAGE_STRING="dynare $VERSION"
make -j"$NTHREADS"
cd "$ROOTDIR"/macOS
cp -L  "$ROOTDIR"/mex/matlab/*                                       "$PKGFILES"/mex/matlab/maci64-9.4-9.7


##
## Make package
##
cd "$ROOTDIR"/macOS/pkg
pkgbuild --root "$PKGFILES" --identifier com.cepremap.dynare --version "$VERSION" --install-location /Applications/Dynare/"$VERSION" "$NAME".pkg
sed "s/VERSION_READ/$VERSION_READ/g" "$ROOTDIR"/macOS/distribution_template.xml > distribution_tmp.xml
sed "s/VERSION_NO_SPACE/$VERSION/g" distribution_tmp.xml > distribution.xml
ln -s "$ROOTDIR"/COPYING "$ROOTDIR"/macOS/
productbuild --distribution distribution.xml --resources "$ROOTDIR"/macOS --package-path ./"$NAME".pkg "$NAME"-new.pkg
rm -f *.xml
rm -rf "$PKGFILES"
rm "$ROOTDIR"/macOS/COPYING
mv "$NAME"-new.pkg "$NAME".pkg
