#!/bin/bash

# Produces Windows packages of Dynare (Windows installer and zip archive).
#
# The binaries are cross compiled for Windows (32/64bits), Octave and MATLAB
# (all supported versions).

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

# Exit on first error and print commands as they are executed
set -ex

# Set root directory
ROOT_DIRECTORY=$(dirname "$(readlink -f "$0")")

# Create TMP folder and make sure it is deleted upon exit
TMP_DIRECTORY=$(mktemp -d)

cleanup()
{
    [[ -z $TMP_DIRECTORY ]] || rm -rf -- "$TMP_DIRECTORY"
}
trap cleanup EXIT

# Set the number of threads
NTHREADS=$(nproc)

# Set Dynare version
if [[ -n $CI_COMMIT_TAG ]]; then
    # Official release tagged through Gitlab
    VERSION=$CI_COMMIT_TAG
    VERSION_SHORT=$VERSION
else
    VERSION=$(grep '^AC_INIT(' ../configure.ac | sed 's/AC_INIT(\[dynare\], \[\(.*\)\])/\1/')
    if [[ -n $CI_COMMIT_SHA ]]; then
        VERSION_SHORT=$VERSION-$CI_COMMIT_SHORT_SHA
        VERSION=$VERSION-$CI_COMMIT_SHA
    elif [[ -d ../.git/ ]]; then
        VERSION_SHORT=$VERSION-$(git rev-parse --short HEAD)
        VERSION=$VERSION-$(git rev-parse HEAD)
    else
        VERSION_SHORT=$VERSION
    fi
fi

BASENAME=dynare-$VERSION_SHORT

# Set directories for dependencies
LIB32="$ROOT_DIRECTORY"/deps/lib32
LIB64="$ROOT_DIRECTORY"/deps/lib64

# Go to source root directory
cd ..

# Autoreconf if needed
[[ -f configure ]] || autoreconf -si

## Compile preprocessor (32-bit), Dynare++ (32-bit) and documentation
./configure --host=i686-w64-mingw32 \
	    --with-boost="$LIB32"/Boost \
	    --with-blas="$LIB32"/OpenBLAS/libopenblas.a \
	    --with-lapack="$LIB32"/OpenBLAS/libopenblas.a \
	    --with-matio="$LIB32"/matIO \
	    --disable-octave \
	    --disable-matlab \
	    PACKAGE_VERSION="$VERSION" \
	    PACKAGE_STRING="dynare $VERSION"
if [[ -z $CI ]]; then
    # If not in Gitlab CI, clean the source and build the doc
    make clean
    make -j"$NTHREADS" pdf html
fi
make -j"$NTHREADS"
i686-w64-mingw32-strip matlab/preprocessor32/dynare_m.exe
i686-w64-mingw32-strip dynare++/src/dynare++.exe

## Compile 64-bit preprocessor
cd preprocessor
make -C src clean # We don't want to clean the doc
./configure --host=x86_64-w64-mingw32 \
	    --with-boost="$LIB64"/Boost \
	    PACKAGE_VERSION="$VERSION" \
	    PACKAGE_STRING="dynare $VERSION"
make -j"$NTHREADS"
x86_64-w64-mingw32-strip src/dynare_m.exe
mkdir -p ../matlab/preprocessor64
mv src/dynare_m.exe ../matlab/preprocessor64

## Define functions for building MEX files

## Note that we do out-of-tree compilation, since we want to do these in
## parallel

# Create Windows 32-bit DLL binaries for MATLAB ≥ R2009b
build_windows_matlab_mex_32 ()
{
    mkdir -p "$TMP_DIRECTORY"/matlab-win32/
    cd "$TMP_DIRECTORY"/matlab-win32/
    "$ROOT_DIRECTORY"/../mex/build/matlab/configure \
                     --host=i686-w64-mingw32 \
		     --with-gsl="$LIB32"/Gsl \
		     --with-matio="$LIB32"/matIO \
		     --with-slicot="$LIB32"/Slicot/without-underscore \
		     --with-matlab="$LIB32"/matlab/R2009b \
		     MATLAB_VERSION=R2009b \
		     MEXEXT=mexw32 \
		     PACKAGE_VERSION="$VERSION" \
		     PACKAGE_STRING="dynare $VERSION"
    make -j"$NTHREADS" all
    i686-w64-mingw32-strip -- **/*.mexw32
    mkdir -p "$ROOT_DIRECTORY"/../mex/matlab/win32-7.9-8.6
    mv -- **/*.mexw32 "$ROOT_DIRECTORY"/../mex/matlab/win32-7.9-8.6
}

# Create Windows 64-bit DLL binaries for MATLAB ≥ R2009b and ≤ R2017b
build_windows_matlab_mex_64_a ()
{
    mkdir -p "$TMP_DIRECTORY"/matlab-win64-a/
    cd "$TMP_DIRECTORY"/matlab-win64-a/
    "$ROOT_DIRECTORY"/../mex/build/matlab/configure \
                     --host=x86_64-w64-mingw32 \
		     --with-gsl="$LIB64"/Gsl \
		     --with-matio="$LIB64"/matIO \
		     --with-slicot="$LIB64"/Slicot/without-underscore \
		     --with-matlab="$LIB64"/matlab/R2009b \
		     MATLAB_VERSION=R2009b \
		     MEXEXT=mexw64 \
		     PACKAGE_VERSION="$VERSION" \
		     PACKAGE_STRING="dynare $VERSION"
    make -j"$NTHREADS" all
    x86_64-w64-mingw32-strip -- **/*.mexw64
    mkdir -p "$ROOT_DIRECTORY"/../mex/matlab/win64-7.9-9.3
    mv -- **/*.mexw64 "$ROOT_DIRECTORY"/../mex/matlab/win64-7.9-9.3
}

# Create Windows 64-bit DLL binaries for MATLAB ≥ R2018a
build_windows_matlab_mex_64_b ()
{
    mkdir -p "$TMP_DIRECTORY"/matlab-win64-b/
    cd "$TMP_DIRECTORY"/matlab-win64-b/
    "$ROOT_DIRECTORY"/../mex/build/matlab/configure \
                     --host=x86_64-w64-mingw32 \
		     --with-gsl="$LIB64"/Gsl \
		     --with-matio="$LIB64"/matIO \
		     --with-slicot="$LIB64"/Slicot/without-underscore \
		     --with-matlab="$LIB64"/matlab/R2018a \
		     MATLAB_VERSION=R2018a \
		     MEXEXT=mexw64 \
		     PACKAGE_VERSION="$VERSION" \
		     PACKAGE_STRING="dynare $VERSION"
    make -j"$NTHREADS" all
    x86_64-w64-mingw32-strip -- **/*.mexw64
    mkdir -p "$ROOT_DIRECTORY"/../mex/matlab/win64-9.4-9.6
    mv -- **/*.mexw64 "$ROOT_DIRECTORY"/../mex/matlab/win64-9.4-9.6
}

# Create Windows DLL binaries for Octave/MinGW (32bit)
build_windows_octave_mex_32 ()
{
    mkdir -p "$TMP_DIRECTORY"/octave-32/
    cd "$TMP_DIRECTORY"/octave-32/
    "$ROOT_DIRECTORY"/../mex/build/octave/configure \
                     --host=i686-w64-mingw32 \
                     --with-gsl="$LIB32"/Gsl \
                     --with-matio="$LIB32"/matIO \
                     --with-slicot="$LIB32"/Slicot/with-underscore \
                     MKOCTFILE="$ROOT_DIRECTORY"/deps/mkoctfile32 \
                     PACKAGE_VERSION="$VERSION" \
                     PACKAGE_STRING="dynare $VERSION"
    make -j"$NTHREADS" all
    i686-w64-mingw32-strip -- **/*.mex
    mkdir -p "$ROOT_DIRECTORY"/../mex/octave/win32
    mv -- **/*.mex "$ROOT_DIRECTORY"/../mex/octave/win32
}

# Create Windows DLL binaries for Octave/MinGW (64bit)
build_windows_octave_mex_64 ()
{
    mkdir -p "$TMP_DIRECTORY"/octave-64/
    cd "$TMP_DIRECTORY"/octave-64/
    "$ROOT_DIRECTORY"/../mex/build/octave/configure \
                     --host=x86_64-w64-mingw32 \
                     --with-gsl="$LIB64"/Gsl \
                     --with-matio="$LIB64"/matIO \
                     --with-slicot="$LIB64"/Slicot/with-underscore \
                     MKOCTFILE="$ROOT_DIRECTORY"/deps/mkoctfile64 \
                     PACKAGE_VERSION="$VERSION" \
                     PACKAGE_STRING="dynare $VERSION"
    make -j"$NTHREADS" all
    x86_64-w64-mingw32-strip -- **/*.mex
    mkdir -p "$ROOT_DIRECTORY"/../mex/octave/win64
    mv -- **/*.mex "$ROOT_DIRECTORY"/../mex/octave/win64
}

## Actually build the MEX files

TASKS=(build_windows_matlab_mex_32 build_windows_matlab_mex_64_a build_windows_matlab_mex_64_b build_windows_octave_mex_32 build_windows_octave_mex_64)
# Reset the number of threads. The mex files for MATLAB/Octave (32-bit and 64-bit) will be built
# in parallel, so we need to account for the number of tasks and lower the value of NTHREADS.
NTHREADS=$((NTHREADS/${#TASKS[@]}))
[[ $NTHREADS -ge 1 ]] || NTHREADS=1 # Ensure that there is at least 1 thread
# Build all the mex files (parallel).
# Some variables and functions need to be available in subshells.
cd "$ROOT_DIRECTORY"
export TMP_DIRECTORY ROOT_DIRECTORY LIB32 LIB64 VERSION NTHREADS
export -f "${TASKS[@]}"
parallel "set -ex;shopt -s globstar;" ::: "${TASKS[@]}"
# Clean up bogus symlinks left by parallel builds of MEX
rm -f ../mex/matlab/*.mexw32 ../mex/matlab/*.mexw64 ../mex/octave/*.mex
# Add supported_octave_version.m (see matlab/dynare.m)
while read -r line
do
    if [[ "$line" =~ OCTAVE_VERSION[[:space:]]*=[[:space:]]*([^[:space:]]+) ]]; then
        OCTAVE_VERSION=${BASH_REMATCH[1]}
        break
    fi
done < "$ROOT_DIRECTORY"/deps/versions.mk
[[ -n $OCTAVE_VERSION ]] || { echo "Can't find OCTAVE_VERSION in versions.mk" >&2; exit 1; }
# shellcheck disable=SC1117
echo -e "function v = supported_octave_version\nv=\"${OCTAVE_VERSION}\";\nend" > ../matlab/supported_octave_version.m

## Create Windows installer
makensis -DVERSION="$VERSION" dynare.nsi
mkdir -p exe
mv dynare-"$VERSION"-win.exe "$ROOT_DIRECTORY"/exe/"$BASENAME"-win.exe

## Create .zip file (for those people that are not allowed to download/execute the installer)

# Set name of the root directory in the ZIP archive
ZIPNAME=dynare-$VERSION
ZIPDIR="$TMP_DIRECTORY"/"$ZIPNAME"
mkdir -p "$ZIPDIR"

cd ..
cp -p NEWS "$ZIPDIR"
cp -p VERSION "$ZIPDIR"
cp -p license.txt "$ZIPDIR"
cp -p windows/README.txt "$ZIPDIR"
cp -pr windows/deps/mingw32 "$ZIPDIR"
cp -pr windows/deps/mingw64 "$ZIPDIR"
mkdir -p "$ZIPDIR"/contrib/ms-sbvar/TZcode
cp -pr contrib/ms-sbvar/TZcode/MatlabFiles "$ZIPDIR"/contrib/ms-sbvar/TZcode
mkdir -p "$ZIPDIR"/contrib/jsonlab
cp -pr contrib/jsonlab/* "$ZIPDIR"/contrib/jsonlab
mkdir "$ZIPDIR"/mex
cp -pr mex/octave/ "$ZIPDIR"/mex
cp -pr mex/matlab/ "$ZIPDIR"/mex
cp -pr matlab "$ZIPDIR"
mkdir -p "$ZIPDIR"/matlab/modules/dseries/externals/x13/windows/32
cp -p windows/deps/lib32/x13as/x13as.exe "$ZIPDIR"/matlab/modules/dseries/externals/x13/windows/32
mkdir -p "$ZIPDIR"/matlab/modules/dseries/externals/x13/windows/64
cp -p windows/deps/lib64/x13as/x13as.exe "$ZIPDIR"/matlab/modules/dseries/externals/x13/windows/64
cp -pr examples "$ZIPDIR"
cp -pr scripts "$ZIPDIR"
mkdir "$ZIPDIR"/dynare++
cp -p dynare++/src/dynare++.exe "$ZIPDIR"/dynare++
mkdir -p "$ZIPDIR"/doc/dynare++
mkdir -p "$ZIPDIR"/doc/dynare-manual.html
cp -pr doc/manual/build/html/* "$ZIPDIR"/doc/dynare-manual.html
cp -p doc/*.pdf "$ZIPDIR"/doc
cp -p doc/manual/build/latex/dynare-manual.pdf "$ZIPDIR"/doc
cp -p preprocessor/doc/macroprocessor/macroprocessor.pdf "$ZIPDIR"/doc
cp -p doc/parallel/parallel.pdf "$ZIPDIR"/doc
cp -p preprocessor/doc/preprocessor/preprocessor.pdf "$ZIPDIR"/doc
cp -p doc/gsa/gsa.pdf "$ZIPDIR"/doc
cp -p dynare++/doc/*.pdf "$ZIPDIR"/doc/dynare++

mkdir -p "$ROOT_DIRECTORY"/zip
cd "$TMP_DIRECTORY"
zip -9 -r "$ROOT_DIRECTORY"/zip/"$BASENAME"-win.zip "$ZIPNAME"
