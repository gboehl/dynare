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

AC_DEFUN([AX_MATLAB_VERSION],
[dnl
AC_REQUIRE([AX_MATLAB])

AC_MSG_CHECKING([for MATLAB version])

if test -n "$MATLAB"; then
  dnl Try to autodetect MATLAB version
  if test -f "${MATLAB}/VersionInfo.xml"; then
     dnl The VersionInfo.xml file is present on all versions since R2017a, on all platforms.
     dnl Extract the version number as x.y, since it is our preferred form, and is
     dnl more robust to future versions.
     dnl NB: brackets in regular expressions are enclosed in other brackets, because this is m4
     MATLAB_VERSION=$(sed -En '/<version>/s/.*>([[0-9]]+\.[[0-9]]+).*/\1/p' "${MATLAB}/VersionInfo.xml")
  elif test -f "${MATLAB}/bin/util/mex/version.txt"; then
     dnl The bin/util/mex/version.txt file is present on Windows and macOS, at least
     dnl since R2009b. It contains the release number (Rnnnnx).
     MATLAB_VERSION=$(cat "${MATLAB}/bin/util/mex/version.txt")
  elif test -f "${MATLAB}/bin/mex" || test -f "${MATLAB}/bin/mexsh"; then
     dnl Works on Linux and macOS until R2018a included. Returns the release number (Rnnnnx).
     dnl Older MATLABs have the version in bin/mex, more recent in bin/mexsh
     dnl NB: brackets in regular expressions are enclosed in other brackets, because this is m4
     MATLAB_VERSION=$(sed -En "/^.*full_ver=/s/^.*full_ver='(R[[^']]+)'.*/\1/p" "${MATLAB}"/bin/mex*)
  fi
fi

if test -n "$MATLAB_VERSION"; then
  dnl Convert a release number (Rnnnnx) into a version number (x.y)
  case $MATLAB_VERSION in
    *2023a | *2023A)
      MATLAB_VERSION="9.14"
      ;;
    *2022b | *2022B)
      MATLAB_VERSION="9.13"
      ;;
    *2022a | *2022A)
      MATLAB_VERSION="9.12"
      ;;
    *2021b | *2021B)
      MATLAB_VERSION="9.11"
      ;;
    *2021a | *2021A)
      MATLAB_VERSION="9.10"
      ;;
    *2020b | *2020B)
      MATLAB_VERSION="9.9"
      ;;
    *2020a | *2020A)
      MATLAB_VERSION="9.8"
      ;;
    *2019b | *2019B)
      MATLAB_VERSION="9.7"
      ;;
    *2019a | *2019A)
      MATLAB_VERSION="9.6"
      ;;
    *2018b | *2018B)
      MATLAB_VERSION="9.5"
      ;;
    *2018a | *2018A)
      MATLAB_VERSION="9.4"
      ;;
    *2017b | *2017B)
      MATLAB_VERSION="9.3"
      ;;
    *2017a | *2017A)
      MATLAB_VERSION="9.2"
      ;;
    *2016b | *2016B)
      MATLAB_VERSION="9.1"
      ;;
    *2016a | *2016A)
      MATLAB_VERSION="9.0"
      ;;
    *2015b | *2015B)
      MATLAB_VERSION="8.6"
      ;;
    *2015a | *2015A)
      MATLAB_VERSION="8.5"
      ;;
    *2014b | *2014B)
      MATLAB_VERSION="8.4"
      ;;
    *2014a | *2014A)
      MATLAB_VERSION="8.3"
      ;;
    *2013b | *2013B)
      MATLAB_VERSION="8.2"
      ;;
    *2013a | *2013A)
      MATLAB_VERSION="8.1"
      ;;
    *2012b | *2012B)
      MATLAB_VERSION="8.0"
      ;;
    *2012a | *2012A)
      MATLAB_VERSION="7.14"
      ;;
    *2011b | *2011B)
      MATLAB_VERSION="7.13"
      ;;
    *2011a | *2011A)
      MATLAB_VERSION="7.12"
      ;;
    *2010b | *2010B)
      MATLAB_VERSION="7.11"
      ;;
    *2010a | *2010A)
      MATLAB_VERSION="7.10"
      ;;
    *2009b | *2009B)
      MATLAB_VERSION="7.9"
      ;;
    *2009a | *2009A)
      MATLAB_VERSION="7.8"
      ;;
    *2008b | *2008B)
      MATLAB_VERSION="7.7"
      ;;
    *2008a | *2008A)
      MATLAB_VERSION="7.6"
      ;;
    *2007b | *2007B)
      MATLAB_VERSION="7.5"
      ;;
    *2007a | *2007A)
      MATLAB_VERSION="7.4"
      ;;
    *2006b | *2006B)
      MATLAB_VERSION="7.3"
      ;;
    *2006a | *2006A)
      MATLAB_VERSION="7.2"
      ;;
    *14SP3 | *14sp3)
      MATLAB_VERSION="7.1"
      ;;
    *14SP2 | *14sp2)
      MATLAB_VERSION="7.0.4"
      ;;
    *14SP1 | *14sp1)
      MATLAB_VERSION="7.0.1"
      ;;
    R14 | r14)
      MATLAB_VERSION="7.0.0"
      ;;
  esac
  dnl Check that we have an x.y version number
  if ! echo "${MATLAB_VERSION}" | grep -qE '^[[0-9.]]+$'; then
      AC_MSG_ERROR([unknown MATLAB version ${MATLAB_VERSION}])
  fi
  AC_MSG_RESULT([${MATLAB_VERSION}])
  ax_matlab_version_ok="yes"
else
  AC_MSG_RESULT([unknown])
  ax_matlab_version_ok="no"
fi
AC_SUBST([MATLAB_VERSION])
])
