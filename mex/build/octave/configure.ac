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
AC_INIT([dynare], [4.7.0])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE([-Wall -Wno-portability foreign])

AX_OCTAVE

test "$ax_enable_octave" != yes && AC_MSG_ERROR([Octave cannot be found])

# Let mkoctfile set the default compilers and flags (except for FC and FCFLAGS which are not supported)
# NB: mkoctfile honors overrides via environment variables
CC=$($MKOCTFILE -p CC)
CXX=$($MKOCTFILE -p CXX)
AR=$($MKOCTFILE -p AR)
RANLIB=$($MKOCTFILE -p RANLIB)
CFLAGS=$($MKOCTFILE -p CFLAGS)
CXXFLAGS=$($MKOCTFILE -p CXXFLAGS)
LDFLAGS="$($MKOCTFILE -p LFLAGS) $($MKOCTFILE -p LDFLAGS)"

AC_CANONICAL_HOST

OCTAVE_VERSION=$($MKOCTFILE -v 2>&1 | sed 's/mkoctfile, version //')
AX_COMPARE_VERSION([$OCTAVE_VERSION], [lt], [5.1.0], [AC_MSG_ERROR([Your Octave is too old, please upgrade to version 5.1.0 at least (or disable Octave support with --disable-octave).])])

AC_PROG_FC
AC_PROG_CC
AC_PROG_CC_C99
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_17
AC_PROG_RANLIB
AX_PROG_LN_S
AC_PROG_MKDIR_P
AM_PROG_AR

AX_CXX11_THREAD

# Check for dlopen(), needed by k_order_perturbation DLL
AC_CHECK_LIB([dl], [dlopen], [LIBADD_DLOPEN="-ldl"], [])
AC_SUBST([LIBADD_DLOPEN])

# Define optional components and their corresponding flags and Automake conditionals
AC_ARG_ENABLE([mex-ms-sbvar], AS_HELP_STRING([--disable-mex-ms-sbvar], [disable compilation of the MS-SBVAR MEX]), [], [enable_mex_ms_sbvar=yes])
AM_CONDITIONAL([ENABLE_MEX_MS_SBVAR], [test "$enable_mex_ms_sbvar" = yes])

AC_ARG_ENABLE([mex-dynare++], AS_HELP_STRING([--disable-mex-dynare++], [disable compilation of MEX based on Dynare++]), [], [enable_mex_dynareplusplus=yes])
AM_CONDITIONAL([ENABLE_MEX_DYNAREPLUSPLUS], [test "$enable_mex_dynareplusplus" = yes])

AC_ARG_ENABLE([mex-kalman-steady-state], AS_HELP_STRING([--disable-mex-kalman-steady-state], [disable compilation of the kalman_steady_state MEX]), [], [enable_mex_kalman_steady_state=yes])
AM_CONDITIONAL([ENABLE_MEX_KALMAN_STEADY_STATE], [test "$enable_mex_kalman_steady_state" = yes])

# Check for GSL, needed by MS-SBVAR
if test "$enable_mex_ms_sbvar" = yes; then
  AX_GSL
  test "$has_gsl" != yes && AC_MSG_ERROR([GSL cannot be found. If you want to skip the compilation of the MS-SBVAR MEX, pass the --disable-mex-ms-sbvar flag.])
fi

# Check for libmatio, needed by MEX files using Dynare++ code, and by ms-sbvar (the latter only under Octave, as an alternative to MATLAB's libmat)
if test "$enable_mex_dynareplusplus" = yes -o "$enable_mex_ms_sbvar" = yes; then
  AX_MATIO
  test "$has_matio" != yes && AC_MSG_ERROR([libmatio cannot be found. If you want to skip the compilation of MS-SBVAR MEX and MEX files based Dynare++, pass the --disable-mex-dynare++ and --disable-mex-ms-sbvar flags.])
fi

# Check for libslicot, needed by kalman_steady_state
if test "$enable_mex_kalman_steady_state" = yes; then
  AX_SLICOT([octave])
  test "$has_slicot" != yes && AC_MSG_ERROR([slicot cannot be found. If you want to skip the compilation of the kalman_steady_state MEX, pass the --disable-mex-kalman-steady-state flag.])
fi

# Check for UMFPACK, needed by bytecode
AC_CHECK_LIB([umfpack], [umfpack_dl_defaults], [LIBADD_UMFPACK="-lumfpack"], [AC_MSG_ERROR([Can't find UMFPACK])])
AC_SUBST([LIBADD_UMFPACK])

# On Windows, we want static linking of the external libraries
case ${host_os} in
  *mingw32*)
    GSL_LIBS="-Wl,-Bstatic $GSL_LIBS -Wl,-Bdynamic"
    LIBADD_MATIO="-Wl,-Bstatic $LIBADD_MATIO -Wl,-Bdynamic"
    LIBADD_SLICOT="-Wl,-Bstatic $LIBADD_SLICOT -Wl,-Bdynamic"
    ;;
esac

# Construct final output message
if test "$enable_mex_dynareplusplus" = yes; then
   BUILD_GENSYLV_KORDER_DYNSIMUL_MEX_OCTAVE="yes"
else
   BUILD_GENSYLV_KORDER_DYNSIMUL_MEX_OCTAVE="no"
fi

if test "$enable_mex_kalman_steady_state" = yes; then
   BUILD_KALMAN_STEADY_STATE_OCTAVE="yes"
else
   BUILD_KALMAN_STEADY_STATE_OCTAVE="no"
fi

if test "$enable_mex_ms_sbvar" = yes; then
   BUILD_MS_SBVAR_MEX_OCTAVE="yes"
else
   BUILD_MS_SBVAR_MEX_OCTAVE="no"
fi

AC_MSG_NOTICE([

Dynare is now configured for building the following components...

Binaries (with "make"):
 MEX files for Octave (except those listed below):                   yes
 Gensylv, k-order and dynare_simul MEX for Octave:                   $BUILD_GENSYLV_KORDER_DYNSIMUL_MEX_OCTAVE
 MS-SBVAR MEX files for Octave:                                      $BUILD_MS_SBVAR_MEX_OCTAVE
 Kalman Steady State MEX file for Octave:                            $BUILD_KALMAN_STEADY_STATE_OCTAVE

])

AC_CONFIG_FILES([Makefile
                 mjdgges/Makefile
                 kronecker/Makefile
                 bytecode/Makefile
                 libdynare++/Makefile
                 gensylv/Makefile
                 libkorder/Makefile
                 k_order_perturbation/Makefile
                 k_order_welfare/Makefile
                 dynare_simul_/Makefile
                 kalman_steady_state/Makefile
                 ms_sbvar/Makefile
                 block_kalman_filter/Makefile
		         sobol/Makefile
		         local_state_space_iterations/Makefile
                 perfect_foresight_problem/Makefile
                 num_procs/Makefile
                 block_trust_region/Makefile
                 disclyap_fast/Makefile])

AC_OUTPUT
