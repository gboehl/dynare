#!/bin/bash
set -exo pipefail

# Creates a dynare-X.Y.mltbx in the current repository, using the settings below.
# Needs to be run from Ubuntu 22.04 LTS, with the needed packages installed.

X13ASVER=1-1-b60
MATLABPATH=/opt/MATLAB/R2024a
# TODO: change size and put white background for better rendering in MATLAB Add-Ons browser
DYNARE_PNG_LOGO=../../preprocessor/doc/logos/dlogo.png

# Prepare temporary workspace and setup cleanup function
tmpdir=$(mktemp -d)
cleanup ()
{
    rm -rf "$tmpdir"
}
trap cleanup EXIT

pushd ../..
meson setup -Dmatlab_path="$MATLABPATH" -Dbuildtype=release -Dprefer_static=true "$tmpdir"/build-matlab-online

cd "$tmpdir"/build-matlab-online
meson compile
meson install --destdir "$tmpdir"
DYNAREVER=$(meson introspect --projectinfo | jq -r '.version')

cd ..
strip usr/local/bin/dynare-preprocessor
strip usr/local/lib/dynare/mex/matlab/*.mexa64

# Get X13 binary from the Census Bureau website
# The binary from Ubuntu has some shared library dependencies, so it is safer to use a static binary
wget -q https://www2.census.gov/software/x-13arima-seats/x13as/unix-linux/program-archives/x13as_ascii-v${X13ASVER}.tar.gz
tar xf x13as_ascii-v${X13ASVER}.tar.gz

# Populate staging area for the zip
cp -pRL usr/local/lib/dynare dynare # -L is needed to dereference the preprocessor symlink
mkdir -p dynare/matlab/dseries/externals/x13/linux/64
cp -p x13as/x13as_ascii dynare/matlab/dseries/externals/x13/linux/64/x13as

# zip dynare
cd dynare
zip -q -r "$tmpdir"/dynare.zip *

# make toolbox
popd
"$MATLABPATH/bin/matlab" -batch "packageDynare('$tmpdir/dynare.zip', '$DYNAREVER', '$DYNARE_PNG_LOGO')"
