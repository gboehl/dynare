#!/bin/bash
set -exo pipefail

# Creates a dynare-X.Y.mltbx in the current repository, using the settings below.
# Needs to be run from Ubuntu 22.04 LTS, with the needed packages installed.

DYNAREVER=5.5
X13ASVER=1-1-b60
LIBLOCATION="/MATLAB Add-Ons/Toolboxes/Dynare/mex/matlab/libs"
MATLABPATH=/opt/MATLAB/R2023b
# MatIO has been recompiled by hand to avoid the dependency on HDF5, which is a nightmare
MATIO_PREFIX=/home/sebastien/usr
# TODO: change size and put white background for better rendering in MATLAB Add-Ons browser
DYNARE_PNG_LOGO=../../preprocessor/doc/logos/dlogo.png

# Prepare temporary workspace and setup cleanup function
tmpdir=$(mktemp -d)
cleanup ()
{
    rm -rf "$tmpdir"
}
trap cleanup EXIT

pushd "$tmpdir"

# Get Dynare
wget -q https://www.dynare.org/release/source/dynare-$DYNAREVER.tar.xz
tar xf dynare-$DYNAREVER.tar.xz

# Build Dynare
cd dynare-$DYNAREVER

# Use static libstdc++ otherwise the one shipped with MATLAB creates problem (since it overrides the systemd-wide one from Ubuntu)
./configure --with-matlab="$MATLABPATH" --with-matio="$MATIO_PREFIX" --disable-octave --disable-doc --disable-dynare++ LDFLAGS=-static-libstdc++
make -j$(nproc)

strip preprocessor/dynare-preprocessor
strip mex/matlab/*.mexa64

# Patch mex files to look into ./lib folder
for f in mex/matlab/*.mexa64
do
    patchelf --set-rpath "$LIBLOCATION" $f
done

# Grab the shared libraries needed
mkdir mex/matlab/libs
for l in libgsl.so.27 libgslcblas.so.0
do
    cp /usr/lib/x86_64-linux-gnu/$l mex/matlab/libs/
    # Patch rpath to find dependencies at runtime
    patchelf --set-rpath "$LIBLOCATION" mex/matlab/libs/$l
done
cp "$MATIO_PREFIX"/lib/libmatio.so.11 mex/matlab/libs/
patchelf --set-rpath "$LIBLOCATION" mex/matlab/libs/libmatio.so.11 # Probably not needed

# Get X13 binary from the Census Bureau website
# The binary from Ubuntu has some shared library dependencies, so it is safer to use a static binary
wget -q https://www2.census.gov/software/x-13arima-seats/x13as/unix-linux/program-archives/x13as_ascii-v${X13ASVER}.tar.gz
tar xf x13as_ascii-v${X13ASVER}.tar.gz
mkdir -p matlab/modules/dseries/externals/x13/linux/64
cp x13as/x13as_ascii matlab/modules/dseries/externals/x13/linux/64/x13as

# zip dynare
zip -q -r "$tmpdir"/dynare.zip contrib/jsonlab contrib/ms-sbvar/TZcode/MatlabFiles examples matlab mex/matlab preprocessor/dynare-preprocessor license.txt NEWS.md README.md VERSION

# make toolbox
popd
"$MATLABPATH/bin/matlab" -batch "packageDynare('$tmpdir/dynare.zip', '$DYNAREVER', '$DYNARE_PNG_LOGO')"
