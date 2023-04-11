dnl Process this file with autoconf to produce a configure script.

dnl Copyright © 2009-2023 Dynare Team
dnl
dnl This file is part of Dynare.
dnl
dnl Dynare is free software: you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.
dnl
dnl Dynare is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

AC_PREREQ([2.62])
AC_INIT([dynare], [6-unstable])
AC_CONFIG_SRCDIR([preprocessor/src/DynareMain.cc])
AM_INIT_AUTOMAKE([1.11 -Wall -Wno-portability foreign no-dist-gzip dist-xz tar-pax])

AC_CONFIG_SUBDIRS([preprocessor])

# For setting EXEEXT, needed by top-level Makefile.am
AC_PROG_CC

AX_PROG_LN_S

AC_PROG_MKDIR_P

# Define optional components and their corresponding flags and Automake conditionals
AC_ARG_ENABLE([doc], AS_HELP_STRING([--disable-doc], [disable compilation of documentation]), [], [enable_doc=yes])
AM_CONDITIONAL([ENABLE_DOC], [test "$enable_doc" = yes])

AC_ARG_ENABLE([matlab], AS_HELP_STRING([--disable-matlab], [disable compilation of MEX files for MATLAB]), [], [enable_matlab=yes])
AM_CONDITIONAL([ENABLE_MATLAB], [test "$enable_matlab" = yes])

AC_ARG_ENABLE([octave], AS_HELP_STRING([--disable-octave], [disable compilation of MEX files for Octave]), [], [enable_octave=yes])
AM_CONDITIONAL([ENABLE_OCTAVE], [test "$enable_octave" = yes])

AC_ARG_ENABLE([org-export], AS_HELP_STRING([--enable-org-export], [enable exporting of Org files (requires Emacs, org-mode and other external programs)]))
AM_CONDITIONAL([ENABLE_ORG_EXPORT], [test -n "$enable_org_export"])

# Check utilities needed for documentation
if test "$enable_doc" = yes; then
  AC_CHECK_PROG([PDFLATEX], [pdflatex], [pdflatex], [no])
  test "$PDFLATEX" = no && AC_MSG_ERROR([pdflatex cannot be found. If you want to skip the compilation of the documentation, pass the --disable-doc flag.])
  AC_CHECK_PROG([BIBTEX], [bibtex], [bibtex], [no])
  test "$BIBTEX" = no && AC_MSG_ERROR([bibtex cannot be found. If you want to skip the compilation of the documentation, pass the --disable-doc flag.])
  AC_CHECK_PROG([SPHINXBUILD], [sphinx-build], [sphinx-build], [no])
  test "$SPHINXBUILD" = no && AC_MSG_ERROR([sphinx-build cannot be found. If you want to skip the compilation of the documentation, pass the --disable-doc flag.])
  AX_LATEX_CLASS([beamer], [ax_latex_have_beamer], [], [AC_MSG_ERROR([beamer cannot be found. If you want to skip the compilation of the documentation, pass the --disable-doc flag.])])
fi

# Check for MATLAB
if test "$enable_matlab" = yes; then
  AC_CONFIG_SUBDIRS([mex/build/matlab])
  AX_MATLAB
  AX_MATLAB_BATCH_OPTIONS
  test "$ax_enable_matlab" != yes && AC_MSG_ERROR([MATLAB cannot be found. If you want to compile Dynare without MATLAB support, pass the --disable-matlab flag.])
fi

# Check for Octave
if test "$enable_octave" = yes; then
  AC_CONFIG_SUBDIRS([mex/build/octave])
  AX_OCTAVE
  test "$ax_enable_octave" != yes && AC_MSG_ERROR([Octave cannot be found. If you want to compile Dynare without Octave support, pass the --disable-octave flag.])
fi

# Construct final output message
if test "$enable_doc" = yes; then
  BUILD_DOC="yes"
else
  BUILD_DOC="no"
fi

if test -n "$enable_org_export"; then
  BUILD_INTERNAL_DOC="yes"
else
  BUILD_INTERNAL_DOC="no"
fi

if test "$enable_matlab" = yes; then
  TESTSUITE_MATLAB="yes"
else
  TESTSUITE_MATLAB="no"
fi

if test "$enable_octave" = yes; then
  TESTSUITE_OCTAVE="yes"
else
  TESTSUITE_OCTAVE="no"
fi

AC_MSG_NOTICE([

Dynare is now configured for building the following components...

Documentation (with "make html pdf"):
 Manual and other documents             $BUILD_DOC
 Dynare internal doc:                   $BUILD_INTERNAL_DOC

Testsuites (run with "make check"):
 Dynare for MATLAB:                     $TESTSUITE_MATLAB
 Dynare for Octave:                     $TESTSUITE_OCTAVE
])

AC_CONFIG_FILES([Makefile
                 VERSION
                 doc/Makefile
                 doc/manual/Makefile
                 doc/manual/utils/version.py
                 doc/parallel/Makefile
                 doc/internals/Makefile
                 doc/gsa/Makefile
                 doc/dseries-and-reporting/Makefile
                 tests/Makefile
                 matlab/dynare_version.m
                 mex/sources/Makefile
])

AC_OUTPUT
