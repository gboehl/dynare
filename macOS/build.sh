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

# Set dependency directory
LIB64="$ROOTDIR"/macOS/deps/lib64

## Hack for statically linking libquadmath, similar to the one used in
## deps/Makefile for several libraries (there is no -static-libquadmath flag,
## see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=46539).
##
## NB: The hack done for Windows does not work for two reasons:
## - the macOS linker is different from GNU ld and does not have the equivalent of -Bstatic/-Bdynamic
## - libgfortran.spec does not include --as-needed on macOS, hence it will link the library anyways
## Also, it does not seem possible to override libgfortran.spec with the --specs option.
GCC_VERSION=$(sed -En "/^c[[:space:]]*=/s/c[[:space:]]*=[[:space:]]*'gcc-([0-9]+)'/\1/p" "$ROOTDIR"/scripts/homebrew-native.ini)
QUADMATH_DIR=$(mktemp -d)
ln -s /usr/local/opt/gcc/lib/gcc/$GCC_VERSION/libquadmath.a $QUADMATH_DIR

##
## Compile Dynare
##
cd "$ROOTDIR"

common_meson_opts=(-Dbuild_for=matlab -Dbuildtype=release -Dprefer_static=true -Dfortran_args="[ '-B', '$LIB64/Slicot/' ]" \
                   -Dc_link_args="[ '-L$QUADMATH_DIR' ]" -Dcpp_link_args="[ '-L$QUADMATH_DIR' ]" -Dfortran_link_args="[ '-L$QUADMATH_DIR' ]" \
                   --native-file scripts/homebrew-native.ini)

# Build for MATLAB ⩾ R2018a
meson setup "${common_meson_opts[@]}" -Dmatlab_path=/Applications/x86_64/MATLAB_R2023b.app build-matlab
meson compile -v -C build-matlab

# Build for MATLAB < R2018a
meson setup "${common_meson_opts[@]}" -Dmatlab_path=/Applications/MATLAB_R2016b.app build-old-matlab
meson compile -v -C build-old-matlab

# If not in CI, build the docs
if [[ -z $CI ]]; then
    meson compile -v -C build-matlab doc
    ln -s build-matlab build-doc
fi

##
## Create package
##

# Determine Dynare version if not passed by an environment variable as in the CI
if [[ -z $VERSION ]]; then
    cd build-matlab
    VERSION=$(meson introspect --projectinfo | sed -En 's/^.*"version": "([^"]*)".*$/\1/p')
    cd ..
fi

# Other useful variables
DATE=$(date +%Y-%m-%d-%H%M)
DATELONG=$(date '+%d %B %Y')

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

NAME=dynare-"$VERSION"
PKGFILES="$ROOTDIR"/macOS/pkg/"$NAME"
mkdir -p \
      "$PKGFILES"/preprocessor \
      "$PKGFILES"/mex/matlab/maci64-8.3-9.3 \
      "$PKGFILES"/mex/matlab/maci64-9.4-23.2 \
      "$PKGFILES"/doc \
      "$PKGFILES"/scripts \
      "$PKGFILES"/contrib/ms-sbvar/TZcode

cp -p  "$ROOTDIR"/NEWS.md                                            "$PKGFILES"
cp -p  "$ROOTDIR"/COPYING                                            "$PKGFILES"
cp -p  "$ROOTDIR"/license.txt                                        "$PKGFILES"

cp -pr "$ROOTDIR"/matlab                                             "$PKGFILES"
cp -p  "$ROOTDIR"/build-matlab/dynare_version.m                      "$PKGFILES"/matlab

cp -pr "$ROOTDIR"/examples                                           "$PKGFILES"

cp -p  "$ROOTDIR"/build-matlab/preprocessor/src/dynare-preprocessor  "$PKGFILES"/preprocessor

# Create backward-compatibility symlink
mkdir -p                                                             "$PKGFILES"/matlab/preprocessor64
ln -sf ../../preprocessor/dynare-preprocessor                        "$PKGFILES"/matlab/preprocessor64/dynare_m

cp -L  "$ROOTDIR"/build-matlab/*.mexmaci64                           "$PKGFILES"/mex/matlab/maci64-9.4-23.2
cp -L  "$ROOTDIR"/build-old-matlab/*.mexmaci64                       "$PKGFILES"/mex/matlab/maci64-8.3-9.3

cp -p  "$ROOTDIR"/scripts/dynare.el                                  "$PKGFILES"/scripts
cp -pr "$ROOTDIR"/contrib/ms-sbvar/TZcode/MatlabFiles                "$PKGFILES"/contrib/ms-sbvar/TZcode
cp -pr "$ROOTDIR"/contrib/jsonlab                                    "$PKGFILES"/contrib

cp     "$ROOTDIR"/build-doc/*.pdf                                    "$PKGFILES"/doc
cp     "$ROOTDIR"/build-doc/preprocessor/doc/*.pdf                   "$PKGFILES"/doc
cp -r  "$ROOTDIR"/build-doc/dynare-manual.html                       "$PKGFILES"/doc

mkdir -p                                                             "$PKGFILES"/matlab/modules/dseries/externals/x13/macOS/64
cp -p  "$ROOTDIR"/macOS/deps/lib64/x13as/x13as                       "$PKGFILES"/matlab/modules/dseries/externals/x13/macOS/64


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
