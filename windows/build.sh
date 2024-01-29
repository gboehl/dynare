#!/bin/bash

# Produces Windows packages of Dynare (executable installer, 7z and zip archives).
#
# The binaries are cross compiled for Windows (64-bit), Octave and MATLAB
# (all supported versions).

# Copyright © 2017-2024 Dynare Team
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

# Exit on first error and print commands as they are executed
set -ex

# Set root directory
ROOT_DIRECTORY=$(dirname "$(readlink -f "$0")")

# Check that build directories do not already exist
[[ -d /tmp/windeps ]] && { echo "Please remove the /tmp/windeps directory" 2>&1; exit 1; }
[[ -d "$ROOT_DIRECTORY"/../build-win-matlab ]] && { echo "Please remove the build-win-matlab directory" 2>&1; exit 1; }
[[ -d "$ROOT_DIRECTORY"/../build-win-octave ]] && { echo "Please remove the build-win-octave directory" 2>&1; exit 1; }

# Create TMP folder and make sure it is deleted upon exit
TMP_DIRECTORY=$(mktemp -d)
cleanup()
{
    [[ -z $TMP_DIRECTORY ]] || rm -rf -- "$TMP_DIRECTORY"
}
trap cleanup EXIT

# Create a directory for dependencies under /tmp.
# Meson does not like when dependencies are under the source tree.
# We use a fixed name to avoid having to regenerate the cross files.
mkdir /tmp/windeps
ln -s "$ROOT_DIRECTORY"/deps/lib64 /tmp/windeps/
ln -s "$ROOT_DIRECTORY"/deps/lib64-msys2 /tmp/windeps/
ln -s "$ROOT_DIRECTORY"/deps/matlab64 /tmp/windeps/
ln -s "$ROOT_DIRECTORY"/deps/mkoctfile64 /tmp/windeps/

# Go to source root directory
cd ..

common_meson_opts=(-Dbuildtype=release --cross-file windows/mingw-cross.ini)

# Create Windows 64-bit DLL binaries for MATLAB ≥ R2018b
meson setup --cross-file windows/mingw-cross-matlab.ini -Dmatlab_path=/tmp/windeps/matlab64/R2018b \
      "${common_meson_opts[@]}" build-win-matlab
meson compile -v -C build-win-matlab

# Create Windows DLL binaries for Octave/MinGW (64bit)
meson setup --cross-file windows/mingw-cross-octave.ini \
      "${common_meson_opts[@]}" build-win-octave
meson compile -v -C build-win-octave

# If not in CI, build the docs
if [[ -z $CI ]]; then
    meson compile -v -C build-win-matlab doc
    ln -sf build-win-matlab build-doc
fi

# Determine Dynare version if not passed by an environment variable as in the CI
if [[ -z $VERSION ]]; then
    cd build-win-matlab
    VERSION=$(meson introspect --projectinfo | sed -En 's/^.*"version": "([^"]*)".*$/\1/p')
    cd ..
fi

# Strip binaries
x86_64-w64-mingw32-strip build-win-matlab/preprocessor/src/dynare-preprocessor.exe
x86_64-w64-mingw32-strip -- build-win-matlab/*.mexw64
x86_64-w64-mingw32-strip -- build-win-octave/*.mex

# Add a preprocessor copy for backward compatibility
mkdir -p matlab/preprocessor64/
cp build-win-matlab/preprocessor/src/dynare-preprocessor.exe matlab/preprocessor64/dynare_m.exe

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
echo -e "function v = supported_octave_version\nv=\"${OCTAVE_VERSION}\";\nend" > matlab/supported_octave_version.m

## Create Windows installer
cd windows
makensis -DVERSION="$VERSION" dynare.nsi
mkdir -p exe
BASENAME=dynare-$VERSION
mv dynare-"$VERSION"-win.exe "$ROOT_DIRECTORY"/exe/"$BASENAME"-win.exe

## Create 7z and zip archives (for people not allowed to download/execute the installer)

# Set name of the root directory in the 7z and zip archives
ZIPDIR="$TMP_DIRECTORY"/"$BASENAME"
mkdir -p "$ZIPDIR"

cd ..
cp -p NEWS.md "$ZIPDIR"
cp -p license.txt "$ZIPDIR"
cp -p windows/README.txt "$ZIPDIR"
cp -pr windows/deps/mingw64 "$ZIPDIR"
mkdir -p "$ZIPDIR"/contrib/ms-sbvar/TZcode
cp -pr contrib/ms-sbvar/TZcode/MatlabFiles "$ZIPDIR"/contrib/ms-sbvar/TZcode
mkdir -p "$ZIPDIR"/mex/matlab/win64-9.5-23.2
cp -p build-win-matlab/*.mexw64 "$ZIPDIR"/mex/matlab/win64-9.5-23.2
mkdir -p "$ZIPDIR"/mex/octave/win64
cp -p build-win-octave/*.mex "$ZIPDIR"/mex/octave/win64
mkdir "$ZIPDIR"/preprocessor
cp -p build-win-matlab/preprocessor/src/dynare-preprocessor.exe "$ZIPDIR"/preprocessor
cp -pr matlab "$ZIPDIR"
cp -p build-win-matlab/dynare_version.m "$ZIPDIR"/matlab
mkdir -p "$ZIPDIR"/matlab/dseries/externals/x13/windows/64
cp -p windows/deps/lib64/x13as/x13as.exe "$ZIPDIR"/matlab/dseries/externals/x13/windows/64
cp -pr examples "$ZIPDIR"
mkdir -p "$ZIPDIR"/scripts
cp -p scripts/dynare.el "$ZIPDIR"/scripts
mkdir -p "$ZIPDIR"/doc
cp -p build-doc/*.pdf "$ZIPDIR"/doc
cp -p build-doc/preprocessor/doc/*.pdf "$ZIPDIR"/doc
cp -pr build-doc/dynare-manual.html "$ZIPDIR"/doc

cd "$TMP_DIRECTORY"

mkdir -p "$ROOT_DIRECTORY"/zip
zip -9 --quiet --recurse-paths "$ROOT_DIRECTORY"/zip/"$BASENAME"-win.zip "$BASENAME"

mkdir -p "$ROOT_DIRECTORY"/7z
7zr a -mx=9 "$ROOT_DIRECTORY"/7z/"$BASENAME"-win.7z "$BASENAME"
