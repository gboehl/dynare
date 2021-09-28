dnl Detect the SLICOT Library.
dnl Called with an argument of either 'matlab' or 'octave', depending
dnl on the configure script from which we're calling it
dnl
dnl AX_SLICOT([matlab])
dnl AX_SLICOT([octave])
dnl
dnl Copyright © 2012-2021 Dynare Team
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

AC_DEFUN([AX_SLICOT],
[
  if test "$1" != matlab && test "$1" != octave; then
    AC_MSG_ERROR([Argument to autoconf slicot macro must be either 'matlab' or 'octave'])
  fi

  AC_ARG_WITH(slicot, AS_HELP_STRING([--with-slicot=DIR], [prefix to SLICOT installation]),
              slicot_prefix="$withval", slicot_prefix="")
  has_slicot=yes

  if test -n "$slicot_prefix"; then
    LDFLAGS_SLICOT="-L$withval/lib"
  else
    LDFLAGS_SLICOT=""
  fi
  my_save_LDFLAGS=$LDFLAGS

  # At this point we should add MATLAB_FCFLAGS to FCFLAGS for Windows (which has -fno-underscoring),
  # but that does not work. The actual underscore test seems to happen at the very beginning of the
  # macro. Hence the modification of FCFLAGS was moved higher (in mex/build/matlab/configure.ac).
  AC_FC_FUNC(sb02od)

  if test "$1" = matlab; then
    LDFLAGS="$LDFLAGS $MATLAB_LDFLAGS_NOMAP $LDFLAGS_SLICOT"

    # Add MATLAB_CFLAGS to get the -fPIC on Linux/x86_64 (otherwise linking fails)
    my_save_CFLAGS=$CFLAGS
    CFLAGS="$CFLAGS $MATLAB_CFLAGS"
    AC_CHECK_LIB([slicot64_pic], [$sb02od], [LIBADD_SLICOT="-lslicot64_pic"], [has_slicot=no], [$MATLAB_LIBS])
    CFLAGS=$my_save_CFLAGS
  else
    LDFLAGS="$LDFLAGS $LDFLAGS_SLICOT"
    # Fallback on libslicot_pic if dynamic libslicot not found
    AC_CHECK_LIB([slicot], [$sb02od], [LIBADD_SLICOT="-lslicot"],
             [
               AC_CHECK_LIB([slicot_pic], [$sb02od], [LIBADD_SLICOT="-lslicot_pic"], [has_slicot=no], [$($MKOCTFILE -p BLAS_LIBS) $($MKOCTFILE -p LAPACK_LIBS)])
             ],
             [$($MKOCTFILE -p BLAS_LIBS) $($MKOCTFILE -p LAPACK_LIBS)])
  fi

  LDFLAGS=$my_save_LDFLAGS
  AC_SUBST(LDFLAGS_SLICOT)
  AC_SUBST(LIBADD_SLICOT)
])
