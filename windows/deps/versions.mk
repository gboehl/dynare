SLICOT_VERSION = 5.0+20101122
X13AS_VERSION = 1.1_B39

OCTAVE_VERSION = 6.2.0
OCTAVE_W32_BUILD =
OCTAVE_W64_BUILD =

MATLAB32_VERSION = 20181204
MATLAB64_VERSION = 20181204


### MSYS2 packages
# Determine the versions by:
# - first running: pacman -Sy
# - and then with appropriate queries using: pacman -Ss <regex>
# Dependencies can be determined using: pacman -Si <pkg>
# File lists can be determined using: pacman -Fl <pkg>
# The same information can be gathered from: https://packages.msys2.org/search

## Build dependencies

# pacman -Ss .*-boost$
MINGW32_BOOST_VERSION = 1.72.0-1
MINGW64_BOOST_VERSION = 1.72.0-1

# pacman -Ss .*-gsl$
MINGW32_GSL_VERSION = 2.6-1
MINGW64_GSL_VERSION = 2.6-1

# pacman -Ss .*-openblas$
MINGW32_OPENBLAS_VERSION = 0.3.7-1
MINGW64_OPENBLAS_VERSION = 0.3.7-1

# pacman -Ss .*-matio$
MINGW32_MATIO_VERSION = 1.5.17-1
MINGW64_MATIO_VERSION = 1.5.17-1

# Dependency of matio (and of the MinGW compiler)
# pacman -Ss .*-zlib$
MINGW32_ZLIB_VERSION = 1.2.11-7
MINGW64_ZLIB_VERSION = 1.2.11-7

# Dependency of matio
# pacman -Ss .*-hdf5$
MINGW32_HDF5_VERSION = 1.10.5-1
MINGW64_HDF5_VERSION = 1.10.5-1

# Dependency of HDF5
# pacman -Ss .*-szip$
MINGW32_SZIP_VERSION = 2.1.1-2
MINGW64_SZIP_VERSION = 2.1.1-2

## MinGW packages for the embedded compiler

# pacman -Ss mingw-w64-.*-gcc$
MINGW32_GCC_VERSION = 9.2.0-2
MINGW64_GCC_VERSION = 9.2.0-2

# pacman -Ss mingw-w64-.*-gmp$
MINGW32_GMP_VERSION = 6.1.2-1
MINGW64_GMP_VERSION = 6.1.2-1

# pacman -Ss mingw-w64-.*-binutils
MINGW32_BINUTILS_VERSION = 2.33.1-1
MINGW64_BINUTILS_VERSION = 2.33.1-1

# pacman -Ss mingw-w64-.*-headers-git
MINGW32_HEADERS_VERSION = 8.0.0.5576.34082b63-1
MINGW64_HEADERS_VERSION = 8.0.0.5576.34082b63-1

# pacman -Ss mingw-w64-.*-crt-git
MINGW32_CRT_VERSION = 8.0.0.5576.34082b63-1
MINGW64_CRT_VERSION = 8.0.0.5576.34082b63-1

# pacman -Ss mingw-w64-.*-winpthreads-git
MINGW32_WINPTHREADS_VERSION = 8.0.0.5574.33e5a2ac-1
MINGW64_WINPTHREADS_VERSION = 8.0.0.5574.33e5a2ac-1

# pacman -Ss mingw-w64-.*-libwinpthread-git
# NB: "thread" is singular here
MINGW32_LIBWINPTHREAD_VERSION = 8.0.0.5574.33e5a2ac-1
MINGW64_LIBWINPTHREAD_VERSION = 8.0.0.5574.33e5a2ac-1
