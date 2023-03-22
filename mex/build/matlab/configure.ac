dnl Process this file with autoconf to produce a configure script.

dnl Copyright © 2009-2022 Dynare Team
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
AC_INIT([dynare], [5.5])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE([-Wall -Wno-portability foreign])

dnl The following must occur before we modify CFLAGS/CXXFLAGS
AC_CANONICAL_HOST

AX_MATLAB
AX_MATLAB_VERSION
AX_DOT_MEXEXT
AX_MATLAB_ARCH
AX_MEXOPTS

if test "$ax_enable_matlab" != yes -o "$ax_matlab_version_ok" != yes -o "$ax_mexopts_ok" != yes; then
  AC_MSG_ERROR([MATLAB cannot be found])
fi

AX_COMPARE_VERSION([$MATLAB_VERSION], [lt], [8.3], [AC_MSG_ERROR([Your MATLAB is too old, please upgrade to version 8.3 (R2014a) at least (or disable MATLAB support with --disable-matlab).])])

case ${host_os} in
  *cygwin*)
    AC_MSG_WARN([You are compiling for the Cygwin target. This means that the MEX files will])
    AC_MSG_WARN([not run from MATLAB, because it cannot load the Cygwin DLL.])
    AC_MSG_WARN([This is probably not what you want. Consider using a MinGW cross-compiler.])
    ;;
esac

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

# Check for GSL, needed by MS-SBVAR MEX
if test "$enable_mex_ms_sbvar" = yes; then
  AX_GSL
  test "$has_gsl" != yes && AC_MSG_ERROR([GSL cannot be found. If you want to skip the compilation of the MS-SBVAR MEX, pass the --disable-mex-ms-sbvar flag.])
fi

# Check for libmatio, needed by MEX files using Dynare++ code
if test "$enable_mex_dynareplusplus" = yes; then
  AX_MATIO
  test "$has_matio" != yes && AC_MSG_ERROR([libmatio cannot be found. If you want to skip the compilation of MEX files based Dynare++, pass the --disable-mex-dynare++ flag.])
fi

# Check for libslicot, needed by kalman_steady_state
if test "$enable_mex_kalman_steady_state" = yes; then
  # FCFLAGS must be temporarily modified, because otherwise -fno-underscoring is not
  # taken into account by AC_FC_FUNC in the AX_SLICOT macro.
  # For some obscure reason, it is necessary to do it at this level and not within the macro.
  my_save_FCFLAGS=$FCFLAGS
  FCFLAGS="$FCFLAGS $MATLAB_FCFLAGS"
  AX_SLICOT([matlab])
  FCFLAGS=$my_save_FCFLAGS
  test "$has_slicot" != yes && AC_MSG_ERROR([slicot cannot be found. If you want to skip the compilation of the kalman_steady_state MEX, pass the --disable-mex-kalman-steady-state flag.])
fi

# On Windows, we want static linking of the external libraries
case ${host_os} in
  *mingw32*)
    GSL_LIBS="-Wl,-Bstatic $GSL_LIBS -Wl,-Bdynamic"
    LIBADD_MATIO="-Wl,-Bstatic $LIBADD_MATIO -Wl,-Bdynamic"
    LIBADD_SLICOT="-Wl,-Bstatic $LIBADD_SLICOT -Wl,-Bdynamic"
    ;;
esac

# Support for M2HTML
AC_ARG_WITH([m2html], AS_HELP_STRING([--with-m2html=DIR], [specify installation directory of M2HTML]), [
M2HTML=$withval
BUILD_M2HTML=yes
], [
M2HTML=
BUILD_M2HTML=no
])
AC_SUBST([M2HTML])
AM_CONDITIONAL([HAVE_M2HTML], [test "x$M2HTML" != "x"])

# Construct final output message
if test "$enable_mex_dynareplusplus" = yes; then
   BUILD_GENSYLV_KORDER_DYNSIMUL_MEX_MATLAB="yes"
else
   BUILD_GENSYLV_KORDER_DYNSIMUL_MEX_MATLAB="no"
fi

if test "$enable_mex_kalman_steady_state" = yes; then
   BUILD_KALMAN_STEADY_STATE_MATLAB="yes"
else
   BUILD_KALMAN_STEADY_STATE_MATLAB="no"
fi

if test "$enable_mex_ms_sbvar" = yes; then
   BUILD_MS_SBVAR_MEX_MATLAB="yes"
else
   BUILD_MS_SBVAR_MEX_MATLAB="no"
fi

AC_MSG_NOTICE([

Dynare is now configured for building the following components...

Binaries (with "make"):
 MEX files for MATLAB (except those listed below):                   yes
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
