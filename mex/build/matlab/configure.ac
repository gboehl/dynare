dnl Process this file with autoconf to produce a configure script.

dnl Copyright © 2009-2019 Dynare Team
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
dnl along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

AC_PREREQ([2.62])
AC_INIT([dynare], [4.6-unstable])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE([-Wall -Wno-portability foreign])

dnl The following must occur before we modify CFLAGS/CXXFLAGS
AC_CANONICAL_HOST

AX_MATLAB

if test "$ax_enable_matlab" = yes; then
  AX_MATLAB_VERSION
  AX_DOT_MEXEXT
  AX_MATLAB_ARCH

  if test "$ax_matlab_version_ok" = yes; then
    AX_MEXOPTS
    AX_COMPARE_VERSION([$MATLAB_VERSION], [lt], [7.9], [AC_MSG_ERROR([Your MATLAB is too old, please upgrade to version 7.9 (R2009b) at least.])])
  fi

  CFLAGS="$MATLAB_CFLAGS"
  CXXFLAGS="$MATLAB_CXXFLAGS"
  FFLAGS="$MATLAB_FFLAGS"
fi

case ${host_os} in
  *cygwin*)
    AC_MSG_WARN([You are compiling for the Cygwin target. This means that the MEX files will])
    AC_MSG_WARN([not run from MATLAB, because it cannot load the Cygwin DLL.])
    AC_MSG_WARN([This is probably not what you want. Consider using a MinGW cross-compiler.])
    ;;
esac

CFLAGS="$CFLAGS -Wall -Wno-parentheses"
FFLAGS="$FFLAGS -Wall"
CXXFLAGS="$CXXFLAGS -Wall -Wno-parentheses"

AC_PROG_F77([gfortran g77 f77])
AC_PROG_CC
AC_PROG_CC_C99 # mjdgges DLL now uses C99 features (variable declared in for loop)
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_14
AC_PROG_RANLIB
AX_PROG_LN_S
AC_PROG_MKDIR_P
AM_PROG_AR

case ${host_os} in
  *darwin*)
    string=$($CXX --version)
    if test "${string#*clang}" != "$string"; then
      CXXFLAGS="$CXXFLAGS -stdlib=libc++"
      MATLAB_LIBS="$MATLAB_LIBS -lc++"
      MATLAB_LDFLAGS="$MATLAB_LDFLAGS -Wl,-syslibroot,$(xcrun -sdk macosx --show-sdk-path)"
    else
      MATLAB_LIBS="$MATLAB_LIBS -lstdc++"
    fi
    ;;
esac

AX_CXX11_THREAD

# Check for dlopen(), needed by k_order_perturbation DLL
AC_CHECK_LIB([dl], [dlopen], [LIBADD_DLOPEN="-ldl"], [])
AC_SUBST([LIBADD_DLOPEN])
# Check for GSL, needed by MS-SBVAR
AX_GSL
AM_CONDITIONAL([HAVE_GSL], [test "$has_gsl" = yes])

# Check for libmatio, needed by MEX files using Dynare++ code
AX_MATIO
AM_CONDITIONAL([HAVE_MATIO], [test "$has_matio" = yes])

# Check for libslicot, needed by kalman_steady_state
AX_SLICOT([matlab])
AM_CONDITIONAL([HAVE_SLICOT], [test "$has_slicot" = yes])

# On Windows, we want static linking of the external libraries
case ${host_os} in
  *mingw32*)
    GSL_LIBS="-Wl,-Bstatic $GSL_LIBS -Wl,-Bdynamic"
    LIBADD_MATIO="-Wl,-Bstatic $LIBADD_MATIO -Wl,-Bdynamic"
    LIBADD_SLICOT="-Wl,-Bstatic $LIBADD_SLICOT -Wl,-Bdynamic"
    ;;
esac

AM_CONDITIONAL([DO_SOMETHING], [test "$ax_enable_matlab" = yes -a "$ax_matlab_version_ok" = yes -a "$ax_mexopts_ok" = yes])

if test "$ax_enable_matlab" = yes -a "$ax_matlab_version_ok" = yes -a "$ax_mexopts_ok" = yes; then
  BUILD_MEX_MATLAB="yes"
else
  BUILD_MEX_MATLAB="no (missing MATLAB, or unknown version, or unknown architecture)"
fi

if test "$ax_enable_matlab" = yes -a "$has_matio" = yes; then
   BUILD_GENSYLV_KORDER_DYNSIMUL_MEX_MATLAB="yes"
else
   BUILD_GENSYLV_KORDER_DYNSIMUL_MEX_MATLAB="no (missing MatIO library)"
fi

if test "$ax_enable_matlab" = yes -a "$ax_matlab_version_ok" = yes -a "$ax_mexopts_ok" = yes -a "$has_slicot" = yes; then
   BUILD_KALMAN_STEADY_STATE_MATLAB="yes"
else
   BUILD_KALMAN_STEADY_STATE_MATLAB="no (missing SLICOT)"
fi

if test "$ax_enable_matlab" = yes -a "$ax_matlab_version_ok" = yes -a "$ax_mexopts_ok" = yes -a "$has_gsl" = yes; then
   BUILD_MS_SBVAR_MEX_MATLAB="yes"
else
   BUILD_MS_SBVAR_MEX_MATLAB="no (missing GSL)"
fi

AC_ARG_ENABLE([openmp], AS_HELP_STRING([--enable-openmp], [use OpenMP for parallelization of some MEX files]), [
  if test "$enable_openmp" = yes; then
    CPPFLAGS="$CPPFLAGS -DUSE_OMP"
    CFLAGS="$CFLAGS -fopenmp"
    CXXFLAGS="$CXXFLAGS -fopenmp"
  fi
])

AC_ARG_WITH([m2html], AS_HELP_STRING([--with-m2html=DIR], [specify installation directory of M2HTML]), [
M2HTML=$withval
BUILD_M2HTML=yes
], [
M2HTML=
BUILD_M2HTML=no
])
AC_SUBST([M2HTML])
AM_CONDITIONAL([HAVE_M2HTML], [test "x$M2HTML" != "x"])

AC_MSG_NOTICE([

Dynare is now configured for building the following components...

Binaries (with "make"):
 MEX files for MATLAB (except those listed below):                   $BUILD_MEX_MATLAB
 Gensylv, k-order and dynare_simul MEX files for MATLAB:             $BUILD_GENSYLV_KORDER_DYNSIMUL_MEX_MATLAB
 MS-SBVAR MEX files for MATLAB:                                      $BUILD_MS_SBVAR_MEX_MATLAB
 Kalman Steady State MEX file for MATLAB:                            $BUILD_KALMAN_STEADY_STATE_MATLAB
 M2HTML documentation:                                               $BUILD_M2HTML

])

AC_CONFIG_FILES([Makefile
                 mjdgges/Makefile
                 kronecker/Makefile
                 bytecode/Makefile
                 libdynare++/Makefile
                 gensylv/Makefile
                 k_order_perturbation/Makefile
                 dynare_simul_/Makefile
                 kalman_steady_state/Makefile
                 ms_sbvar/Makefile
                 block_kalman_filter/Makefile
	         sobol/Makefile
		 local_state_space_iterations/Makefile])

AC_OUTPUT
