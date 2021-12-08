dnl Process this file with autoconf to produce a configure script.

dnl Copyright © 2009-2021 Dynare Team
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

AC_PROG_CC
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_17

AC_CANONICAL_HOST
case ${host_os} in
  *mingw32*)
    # On Windows, we want Dynare++ to be statically linked
    AM_LDFLAGS="-static"
    AC_SUBST([AM_LDFLAGS])
    ;;
  *cygwin*)
    AC_MSG_WARN([You are compiling for the Cygwin target. This means that the preprocessor will])
    AC_MSG_WARN([not run from MATLAB unless you add the Cygwin DLL to the path.])
    AC_MSG_WARN([This is probably not what you want. Consider using a MinGW cross-compiler.])
    if test -z "$F77"; then
      # On Cygwin 1.7, g77 comes has version 3, and is not compatible with default gcc/g++ which has version 4
      # And by default, the AC_PROG_F77 will pick up g77 if it is present (even if gfortran is also here)
      F77=gfortran
    fi
    ;;
esac

# Use C++ for testing headers
AC_LANG([C++])

AM_CXXFLAGS="-Wall -Wno-parentheses -Wold-style-cast"
AC_SUBST([AM_CXXFLAGS])

# If default 'ar' is not available, try to find one with a host prefix (see ticket #145)
AC_CHECK_PROGS([AR], [ar ${host_alias}-ar])

AC_PROG_RANLIB
AM_PROG_AR

AX_PROG_LN_S

AC_PROG_MKDIR_P

AX_CXX11_THREAD

AM_PROG_LEX

AC_CHECK_PROG([YACC], [bison], [bison])
if test -z "$YACC"; then
  unset YACC # AM_MISSING_PROG needs an unset variable: an empty variable won't do
  AM_MISSING_PROG([YACC], [bison])
fi

# Define optional components and their corresponding flags and Automake conditionals
AC_ARG_ENABLE([doc], AS_HELP_STRING([--disable-doc], [disable compilation of documentation]), [], [enable_doc=yes])
AM_CONDITIONAL([ENABLE_DOC], [test "$enable_doc" = yes])

AC_ARG_ENABLE([dynare++], AS_HELP_STRING([--disable-dynare++], [disable compilation of Dynare++]), [], [enable_dynareplusplus=yes])
AM_CONDITIONAL([ENABLE_DYNAREPLUSPLUS], [test "$enable_dynareplusplus" = yes])

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

# Check for BLAS, LAPACK and MATIO, needed by Dynare++
if test "$enable_dynareplusplus" = yes; then
  AC_PROG_F77
  AC_F77_LIBRARY_LDFLAGS
  case ${host_os} in
    *darwin*)
      # * OS X doesn't need to link against the Fortran libraries because
      #   everything is contained within the vecLib framework.
      FLIBS=''
      ;;
  esac
  AX_BLAS
  AX_LAPACK
  AX_MATIO
  if test "$ax_blas_ok" != yes -o "$ax_lapack_ok" != yes -o "$has_matio" != yes; then
    AC_MSG_ERROR([Some dependencies of Dynare++ cannot be found. If you want to skip the compilation of Dynare++, pass the --disable-dynare++ flag.])
  fi
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
if test "$enable_dynareplusplus" = yes; then
  BUILD_DYNAREPLUSPLUS="yes"
  TESTSUITE_DYNAREPLUSPLUS="yes"
else
  BUILD_DYNAREPLUSPLUS="no"
  TESTSUITE_DYNAREPLUSPLUS="no"
fi

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

Binaries (with "make"):
 Dynare++:                              $BUILD_DYNAREPLUSPLUS

Documentation (with "make html pdf"):
 Manual and other documents             $BUILD_DOC
 Dynare internal doc:                   $BUILD_INTERNAL_DOC

Testsuites (run with "make check"):
 Dynare for MATLAB:                     $TESTSUITE_MATLAB
 Dynare for Octave:                     $TESTSUITE_OCTAVE
 Dynare++:                              $TESTSUITE_DYNAREPLUSPLUS
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
                 dynare++/Makefile
                 dynare++/parser/cc/Makefile
                 dynare++/sylv/Makefile
                 dynare++/sylv/cc/Makefile
                 dynare++/sylv/testing/Makefile
                 dynare++/tl/Makefile
                 dynare++/tl/cc/Makefile
                 dynare++/tl/testing/Makefile
                 dynare++/doc/Makefile
                 dynare++/utils/cc/Makefile
                 dynare++/integ/Makefile
                 dynare++/integ/cc/Makefile
                 dynare++/integ/src/Makefile
                 dynare++/integ/testing/Makefile
                 dynare++/kord/Makefile
                 dynare++/src/Makefile
                 dynare++/tests/Makefile
                 mex/sources/Makefile
])

AC_OUTPUT
