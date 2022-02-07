#!/bin/bash

# On a Debian system, install the packages needed for Windows
# cross-compilation, and also setup the cross-compiler alternatives.

# Copyright © 2017-2022 Dynare Team
#
# This file is part of Dynare.
#
# Dynare is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Dynare is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

[[ $(id -u) == 0 ]] || { echo "You must be root" >&2; exit 1; }

PACKAGES=(make xz-utils p7zip bzip2 zip zstd patch wget autoconf automake
          libtool mingw-w64 gfortran-mingw-w64 parallel flex libfl-dev bison texlive
          texlive-publishers texlive-latex-extra texlive-science
          texlive-fonts-extra lmodern python3-sphinx latexmk nsis)

apt install "${PACKAGES[@]}"

# Configure MinGW to use the POSIX threading model (needed for C++11 threads in
# Dynare++, see /usr/share/doc/gcc-mingw-w64-base/README.Debian)
update-alternatives --set x86_64-w64-mingw32-gfortran /usr/bin/x86_64-w64-mingw32-gfortran-posix
update-alternatives --set x86_64-w64-mingw32-gcc /usr/bin/x86_64-w64-mingw32-gcc-posix
update-alternatives --set x86_64-w64-mingw32-g++ /usr/bin/x86_64-w64-mingw32-g++-posix
