#!/bin/bash

# Produces Windows packages of Dynare (executable installer, 7z and zip archives).
#
# The binaries are cross compiled for Windows (32/64bits), Octave and MATLAB
# (all supported versions).

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

# Set Dynare version, if not already set by Gitlab CI
if [[ -z $VERSION ]]; then
    VERSION=$(grep '^AC_INIT(' ../configure.ac | sed 's/AC_INIT(\[dynare\], \[\(.*\)\])/\1/')
    if [[ -d ../.git/ ]]; then
        VERSION=$VERSION-$(git rev-parse --short HEAD)
    fi
fi

BASENAME=dynare-$VERSION

# Set directories for dependencies
LIB32="$ROOT_DIRECTORY"/deps/lib32
LIB64="$ROOT_DIRECTORY"/deps/lib64
LIB32_MSYS2="$ROOT_DIRECTORY"/deps/lib32-msys2
LIB64_MSYS2="$ROOT_DIRECTORY"/deps/lib64-msys2

# Go to source root directory
cd ..

# Autoreconf if needed
[[ -f configure ]] || autoreconf -si

## Compile preprocessor (64-bit), Dynare++ (64-bit) and documentation
./configure --host=x86_64-w64-mingw32 \
	    --with-boost="$LIB64_MSYS2" \
	    --with-blas="$LIB64_MSYS2"/lib/libopenblas.a \
	    --with-lapack="$LIB64_MSYS2"/lib/libopenblas.a \
	    --with-matio="$LIB64_MSYS2" \
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
x86_64-w64-mingw32-strip matlab/preprocessor64/dynare_m.exe
x86_64-w64-mingw32-strip dynare++/src/dynare++.exe
mkdir -p dynare++/64-bit/
cp dynare++/src/dynare++.exe dynare++/64-bit/

## Compile 32-bit preprocessor and Dynare++
# We do not want to clean the doc
for d in preprocessor/src dynare++/integ dynare++/kord dynare++/parser/cc dynare++/src dynare++/sylv dynare++/tl dynare++/utils/cc; do
    make -C "$d" clean
done
./configure --host=i686-w64-mingw32 \
	    --with-boost="$LIB32_MSYS2" \
	    --with-blas="$LIB32_MSYS2"/lib/libopenblas.a \
	    --with-lapack="$LIB32_MSYS2"/lib/libopenblas.a \
	    --with-matio="$LIB32_MSYS2" \
	    --disable-octave \
	    --disable-matlab \
	    PACKAGE_VERSION="$VERSION" \
	    PACKAGE_STRING="dynare $VERSION"
make -j"$NTHREADS"
i686-w64-mingw32-strip matlab/preprocessor32/dynare_m.exe
i686-w64-mingw32-strip dynare++/src/dynare++.exe
mkdir -p dynare++/32-bit/
cp dynare++/src/dynare++.exe dynare++/32-bit/

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
		     --with-gsl="$LIB32_MSYS2" \
		     --with-matio="$LIB32_MSYS2" \
		     --with-slicot="$LIB32"/Slicot/without-underscore \
		     --with-matlab="$ROOT_DIRECTORY"/deps/matlab32/R2009b \
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
		     --with-gsl="$LIB64_MSYS2" \
		     --with-matio="$LIB64_MSYS2" \
		     --with-slicot="$LIB64"/Slicot/without-underscore \
		     --with-matlab="$ROOT_DIRECTORY"/deps/matlab64/R2009b \
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
		     --with-gsl="$LIB64_MSYS2" \
		     --with-matio="$LIB64_MSYS2" \
		     --with-slicot="$LIB64"/Slicot/without-underscore \
		     --with-matlab="$ROOT_DIRECTORY"/deps/matlab64/R2018a \
		     MATLAB_VERSION=R2018a \
		     MEXEXT=mexw64 \
		     PACKAGE_VERSION="$VERSION" \
		     PACKAGE_STRING="dynare $VERSION"
    make -j"$NTHREADS" all
    x86_64-w64-mingw32-strip -- **/*.mexw64
    mkdir -p "$ROOT_DIRECTORY"/../mex/matlab/win64-9.4-9.7
    mv -- **/*.mexw64 "$ROOT_DIRECTORY"/../mex/matlab/win64-9.4-9.7
}

# Create Windows DLL binaries for Octave/MinGW (32bit)
build_windows_octave_mex_32 ()
{
    mkdir -p "$TMP_DIRECTORY"/octave-32/
    cd "$TMP_DIRECTORY"/octave-32/
    "$ROOT_DIRECTORY"/../mex/build/octave/configure \
                     --host=i686-w64-mingw32 \
                     --with-gsl="$LIB32_MSYS2" \
                     --with-matio="$LIB32_MSYS2" \
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
                     --with-gsl="$LIB64_MSYS2" \
                     --with-matio="$LIB64_MSYS2" \
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
export TMP_DIRECTORY ROOT_DIRECTORY LIB32 LIB32_MSYS2 LIB64 LIB64_MSYS2 VERSION NTHREADS
export -f "${TASKS[@]}"
parallel "set -ex;shopt -s globstar;" ::: "${TASKS[@]}"

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

## Create 7z and zip archives (for people not allowed to download/execute the installer)

# Set name of the root directory in the 7z and zip archives
ZIPNAME=dynare-$VERSION
ZIPDIR="$TMP_DIRECTORY"/"$ZIPNAME"
mkdir -p "$ZIPDIR"

cd ..
cp -p NEWS.md "$ZIPDIR"
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
mkdir -p "$ZIPDIR"/scripts
cp -p scripts/dynare.el "$ZIPDIR"/scripts
mkdir "$ZIPDIR"/dynare++
cp -pr dynare++/32-bit/ dynare++/64-bit/ dynare++/dynare_simul/dynare_simul.m "$ZIPDIR"/dynare++
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

cd "$TMP_DIRECTORY"

mkdir -p "$ROOT_DIRECTORY"/zip
zip -9 --quiet --recurse-paths "$ROOT_DIRECTORY"/zip/"$BASENAME"-win.zip "$ZIPNAME"

mkdir -p "$ROOT_DIRECTORY"/7z
7zr a -mx=9 "$ROOT_DIRECTORY"/7z/"$BASENAME"-win.7z "$ZIPNAME"
