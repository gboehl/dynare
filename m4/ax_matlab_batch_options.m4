dnl Copyright © 2019 Dynare Team
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

AC_DEFUN([AX_MATLAB_BATCH_OPTIONS],
[dnl
AC_REQUIRE([AX_MATLAB])
AC_REQUIRE([AX_MATLAB_ARCH])
AC_REQUIRE([AX_MATLAB_VERSION])

AX_COMPARE_VERSION([$MATLAB_VERSION], [ge], [9.6],
  [
    if test "${MATLAB_ARCH}" = win32 -o "${MATLAB_ARCH}" = win64; then
      MATLAB_BATCH_OPTIONS='-noFigureWindows -batch'
    else
      MATLAB_BATCH_OPTIONS='-nodisplay -batch'
    fi
  ],
  [
    if test "${MATLAB_ARCH}" = win32 -o "${MATLAB_ARCH}" = win64; then
      MATLAB_BATCH_OPTIONS='-nosplash -automation -wait -sd "$(CURDIR)" -r'
    else
      MATLAB_BATCH_OPTIONS='-nosplash -nodisplay -r'
    fi
  ])

AC_SUBST([MATLAB_BATCH_OPTIONS])
])
