#!/bin/bash
set -ex

VERSION=4.4.3
TOP_DIR=/Users/Houtan/Documents/DYNARE/PACKAGES
TOP_DYN_DIR=$TOP_DIR/dynare-$VERSION

INSTALLDIRNAME=dynare-$VERSION-osx
rm -rf $INSTALLDIRNAME
mkdir $INSTALLDIRNAME
INSTALLDIR=$TOP_DIR/$INSTALLDIRNAME

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
mkdir -p "$INSTALLDIR/mex/octave"
mkdir -p "$INSTALLDIR/mex/matlab/osx64"
mkdir "$INSTALLDIR/mex/matlab/osx32-7.4"
mkdir "$INSTALLDIR/mex/matlab/osx32-7.5-7.11"

# top level
cp $TOP_DYN_DIR/scripts/dynare.el                                $INSTALLDIR/scripts
cp $TOP_DYN_DIR/license.txt                                      $INSTALLDIR
cp $TOP_DYN_DIR/NEWS                                             $INSTALLDIR

# TZ Matlab
cp -r $TOP_DYN_DIR/contrib/ms-sbvar/TZcode/MatlabFiles           $INSTALLDIR/contrib/ms-sbvar/TZcode

# examples
cp -r $TOP_DYN_DIR/examples                                      $INSTALLDIR

##########################################################
# FIRST BUILD 32 BIT EVERYTHING, 32 BIT MATLAB < 7.5 MEX #
##########################################################
./configure FFLAGS='-isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.6.sdk -mmacosx-version-min=10.6 -arch i386' CPPFLAGS='-isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.6.sdk -mmacosx-version-min=10.6 -arch i386' LDFLAGS='-isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.6.sdk -mmacosx-version-min=10.6 -arch i386' --disable-octave --with-matlab=/Applications/MATLAB_OLD/R2007a MATLAB_VERSION=7.4 --with-gsl=/usr/local32 --with-slicot=/usr/local32 --with-matio=/usr/local32
cd $TOP_DYN_DIR/doc
texi2dvi --pdf --batch --build-dir=dynare.t2p dynare.texi

cd $TOP_DYN_DIR
make pdf

cd $TOP_DYN_DIR/preprocessor
make

cd $TOP_DYN_DIR/dynare++
make

cd $TOP_DYN_DIR/mex/build/matlab
make

# Matlab
# Must come after configure because matlab/dynare_version.m is created by configure script
cp -r $TOP_DYN_DIR/matlab                                        $INSTALLDIR

########################
# MAKE BULK OF PACKAGE #
########################
# compiled preprocessor
cp $TOP_DYN_DIR/preprocessor/dynare_m                            $INSTALLDIR/matlab

# Matlab Mex
cp $TOP_DYN_DIR/mex/matlab/*                                     $INSTALLDIR/mex/matlab/osx32-7.4

# dynare++
cp $TOP_DYN_DIR/dynare++/src/dynare++                                               $INSTALLDIR/dynare++
cp $TOP_DYN_DIR/dynare++/extern/matlab/dynare_simul.m                               $INSTALLDIR/dynare++

# doc
cp $TOP_DYN_DIR/doc/bvar-a-la-sims.pdf                                              $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/dr.pdf                                                          $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/dynare.pdf                                                      $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/guide.pdf                                                       $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/macroprocessor/macroprocessor.pdf                               $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/parallel/parallel.pdf                                           $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/preprocessor/preprocessor.pdf                                   $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/userguide/UserGuide.pdf                                         $INSTALLDIR/doc
cp $TOP_DYN_DIR/doc/gsa/gsa.pdf                                                     $INSTALLDIR/doc

# doc (dynare++)
cp $TOP_DYN_DIR/dynare++/doc/dynare++-tutorial.pdf                                  $INSTALLDIR/doc/dynare++
cp $TOP_DYN_DIR/dynare++/doc/dynare++-ramsey.pdf                                    $INSTALLDIR/doc/dynare++
cp $TOP_DYN_DIR/dynare++/sylv/sylvester.pdf                                         $INSTALLDIR/doc/dynare++
cp $TOP_DYN_DIR/dynare++/tl/cc/tl.pdf                                               $INSTALLDIR/doc/dynare++
cp $TOP_DYN_DIR/dynare++/integ/cc/integ.pdf                                         $INSTALLDIR/doc/dynare++
cp $TOP_DYN_DIR/dynare++/kord/kord.pdf                                              $INSTALLDIR/doc/dynare++


##############################################
# RETURN TO BUILD 32 BIT MATLAB 7.5 & UP MEX #
##############################################
cd $TOP_DYN_DIR/mex/build/matlab
make clean
./configure FFLAGS='-isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.6.sdk -mmacosx-version-min=10.6 -arch i386' CPPFLAGS='-isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.6.sdk -mmacosx-version-min=10.6 -arch i386' LDFLAGS='-isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.6.sdk -mmacosx-version-min=10.6 -arch i386' --with-matlab=/Applications/MATLAB_OLD/MATLAB_R2009b_32bit/MATLAB_R2009b.app MATLAB_VERSION=7.9 MEXEXT='mexmaci' --with-slicot=/usr/local32 --with-matio=/usr/local32 --with-gsl=/usr/local32
make

# Matlab Mex
cp $TOP_DYN_DIR/mex/matlab/*                                                         $INSTALLDIR/mex/matlab/osx32-7.5-7.11

#####################################
# RETURN TO BUILD 64 BIT MATLAB MEX #
#####################################
cd $TOP_DYN_DIR/mex/build/matlab
make clean
./configure --with-matlab=/Applications/MATLAB_OLD/MATLAB_R2009b.app MATLAB_VERSION=7.9 MEXEXT=mexmaci64 --with-matio=/usr/localStatic --with-gsl=/usr/localStatic
make

# Matlab Mex
cp $TOP_DYN_DIR/mex/matlab/*                                                           $INSTALLDIR/mex/matlab/osx64

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
