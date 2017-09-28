#!/bin/bash
set -ex

VERSION=4.5.0
TOP_DIR=/Users/houtanb/Documents/DYNARE/package
TOP_DYN_DIR=$TOP_DIR/dynare-$VERSION

INSTALLDIRNAME=dynare-$VERSION-osx
INSTALLDIR=$TOP_DIR/$INSTALLDIRNAME
rm -rf $INSTALLDIR
mkdir $INSTALLDIR

########################
# UPDATE DYNARE SOURCE #
########################
cd $TOP_DYN_DIR

########################
# BEGIN MAKING PACKAGE #
########################
# create directories
mkdir "$INSTALLDIR/scripts"
mkdir "$INSTALLDIR/dynare++"
mkdir -p "$INSTALLDIR/doc/dynare++"
mkdir "$INSTALLDIR/doc/dynare.html"
mkdir -p "$INSTALLDIR/contrib/ms-sbvar/TZcode"
mkdir -p "$INSTALLDIR/mex/matlab/osx"

# top level
cp $TOP_DYN_DIR/scripts/dynare.el                                $INSTALLDIR/scripts
cp $TOP_DYN_DIR/NEWS                                             $INSTALLDIR
cp $TOP_DYN_DIR/COPYING                                          $INSTALLDIR
cp $TOP_DYN_DIR/license.txt                                      $INSTALLDIR

# TZ Matlab
cp -r $TOP_DYN_DIR/contrib/ms-sbvar/TZcode/MatlabFiles           $INSTALLDIR/contrib/ms-sbvar/TZcode

# examples
cp -r $TOP_DYN_DIR/examples                                      $INSTALLDIR

#############
# CONFIGURE #
#############
./configure FFLAGS='-mmacosx-version-min=10.8' CPPFLAGS='-mmacosx-version-min=10.8' LDFLAGS='-mmacosx-version-min=10.8' \
            --with-matlab=/Volumes/Storage/MATLAB/MATLAB_R2009b.app MATLAB_VERSION=7.9 MEXEXT=mexmaci64 --disable-octave \
            --with-matio=/usr/local/static --with-gsl=/usr/local/static --with-slicot=/usr/local/static
make clean

###########
# COMPILE #
###########
cd $TOP_DYN_DIR/doc
texi2dvi --pdf --batch --build-dir=dynare.t2p dynare.texi

cd $TOP_DYN_DIR
make
make pdf

###########################
# COMPOSE BULK OF PACKAGE #
###########################

# Matlab
# Must come after configure because matlab/dynare_version.m is created by configure script
cp -r $TOP_DYN_DIR/matlab                                        $INSTALLDIR

# compiled preprocessor
cp $TOP_DYN_DIR/preprocessor/dynare_m                            $INSTALLDIR/matlab

# dynare++
cp $TOP_DYN_DIR/dynare++/src/dynare++                            $INSTALLDIR/dynare++
cp $TOP_DYN_DIR/dynare++/extern/matlab/dynare_simul.m            $INSTALLDIR/dynare++

# doc
cp $TOP_DYN_DIR/doc/bvar-a-la-sims.pdf                           $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/dr.pdf                                       $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/dynare.pdf                                   $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/guide.pdf                                    $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/macroprocessor/macroprocessor.pdf            $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/parallel/parallel.pdf                        $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/preprocessor/preprocessor.pdf                $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/gsa/gsa.pdf                                  $INSTALLDIR/doc

# doc (dynare++)
cp $TOP_DYN_DIR/dynare++/doc/dynare++-tutorial.pdf               $INSTALLDIR/doc/dynare++
cp $TOP_DYN_DIR/dynare++/doc/dynare++-ramsey.pdf                 $INSTALLDIR/doc/dynare++
cp $TOP_DYN_DIR/dynare++/sylv/sylvester.pdf                      $INSTALLDIR/doc/dynare++
cp $TOP_DYN_DIR/dynare++/tl/cc/tl.pdf                            $INSTALLDIR/doc/dynare++
cp $TOP_DYN_DIR/dynare++/integ/cc/integ.pdf                      $INSTALLDIR/doc/dynare++
cp $TOP_DYN_DIR/dynare++/kord/kord.pdf                           $INSTALLDIR/doc/dynare++

# Matlab Mex
cp $TOP_DYN_DIR/mex/matlab/*                                     $INSTALLDIR/mex/matlab/osx

# clean everything
cd $TOP_DYN_DIR

# remove .DS_Store files
cd $INSTALLDIR
find . -name *.DS_Store -type f -exec rm {} \;

# Change permissions
chmod -R g+w $INSTALLDIR

echo DONE
# NEED TO BUILD DYNARE.HTML DOCUMENTION ON DEBIAN
# AND INCLUDE IN DISTRIBUTION BY HAND
