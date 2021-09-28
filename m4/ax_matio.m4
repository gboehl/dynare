dnl Detect the MATIO Library.
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

AC_DEFUN([AX_MATIO],
[
AC_ARG_WITH(matio, AS_HELP_STRING([--with-matio=DIR], [prefix to MATIO installation]),
            matio_prefix="$withval", matio_prefix="")

  has_matio=yes

  if test -n "$matio_prefix"; then
    CPPFLAGS_MATIO="-I$withval/include"
    LDFLAGS_MATIO="-L$withval/lib"
  else
    CPPFLAGS_MATIO=""
    LDFLAGS_MATIO=""
  fi

  ac_save_CPPFLAGS="$CPPFLAGS"
  ac_save_LDFLAGS="$LDFLAGS"
  ac_save_LIBS="$LIBS"

  LIBADD_MATIO=""
  CPPFLAGS="$CPPFLAGS_MATIO $CPPFLAGS"
  LDFLAGS="$LDFLAGS_MATIO $LDFLAGS"

  dnl Workaround for the matio from RHEL 6 + EPEL 6
  dnl If detected, libz and libhdf5 are added to LIBS, used for matio test
  LIBS=""
  AC_CHECK_LIB([z], [compress])
  dnl szip is needed under MSYS2
  AC_CHECK_LIB([szip], [SZ_Compress])
  dnl szip is needed for static linking on macOS (it's called libsz on macOS)
  AC_CHECK_LIB([sz], [SZ_Compress])
  AC_CHECK_LIB([hdf5], [H5Fcreate])

  AC_CHECK_HEADER([matio.h], [], [has_matio=no])
  AC_CHECK_LIB([matio], [Mat_Open], [LIBADD_MATIO="-lmatio $LIBS"], [has_matio=no])

  CPPFLAGS="$ac_save_CPPFLAGS"
  LDFLAGS="$ac_save_LDFLAGS"
  LIBS="$ac_save_LIBS"

  AC_SUBST(CPPFLAGS_MATIO)
  AC_SUBST(LIBADD_MATIO)
  AC_SUBST(LDFLAGS_MATIO)
])
