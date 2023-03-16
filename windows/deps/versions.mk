SLICOT_VERSION = 5.0+20101122
X13AS_VERSION = 1-1-b59

OCTAVE_VERSION = 8.1.0
OCTAVE_W64_BUILD =

MATLAB64_VERSION = 20200930


### MSYS2 packages
# Determine the versions by:
# - first running: pacman -Sy
# - and then with appropriate queries using: pacman -Ss <regex>
# Dependencies can be determined using: pacman -Si <pkg>
# File lists can be determined using: pacman -Fl <pkg>
# The same information can be gathered from: https://packages.msys2.org/search

## Build dependencies

# pacman -Ss mingw-w64-x86_64-boost
MINGW64_BOOST_VERSION = 1.81.0-6

# pacman -Ss mingw-w64-x86_64-gsl
MINGW64_GSL_VERSION = 2.7.1-1

# pacman -Ss mingw-w64-x86_64-matio
MINGW64_MATIO_VERSION = 1.5.23-4

# Dependency of matio (and of the MinGW compiler)
# pacman -Ss mingw-w64-x86_64-zlib
MINGW64_ZLIB_VERSION = 1.2.13-3

# Dependency of matio
# pacman -Ss mingw-w64-x86_64-hdf5
MINGW64_HDF5_VERSION = 1.12.2-2

# Dependency of HDF5 (provides szip library)
# pacman -Ss mingw-w64-x86_64-libaec
MINGW64_LIBAEC_VERSION = 1.0.6-2

## MinGW packages for the embedded compiler

# pacman -Ss mingw-w64-x86_64-gcc$
MINGW64_GCC_VERSION = 12.2.0-10

# pacman -Ss mingw-w64-x86_64-gmp
MINGW64_GMP_VERSION = 6.2.1-5

# pacman -Ss mingw-w64-x86_64-binutils
MINGW64_BINUTILS_VERSION = 2.40-2

# pacman -Ss mingw-w64-x86_64-headers-git
MINGW64_HEADERS_VERSION = 10.0.0.r234.g283e5b23a-1

# pacman -Ss mingw-w64-x86_64-crt-git
MINGW64_CRT_VERSION = 10.0.0.r234.g283e5b23a-1

# pacman -Ss mingw-w64-x86_64-winpthreads-git
MINGW64_WINPTHREADS_VERSION = 10.0.0.r234.g283e5b23a-1

# pacman -Ss mingw-w64-x86_64-zstd
MINGW64_ZSTD_VERSION = 1.5.4-1

# pacman -Ss mingw-w64-x86_64-isl
MINGW64_ISL_VERSION = 0.25-1

# pacman -Ss mingw-w64-x86_64-mpc
MINGW64_MPC_VERSION = 1.3.1-1

# pacman -Ss mingw-w64-x86_64-mpfr
MINGW64_MPFR_VERSION = 4.2.0-1

# pacman -Ss mingw-w64-x86_64-libiconv
MINGW64_LIBICONV_VERSION = 1.17-3

# pacman -Ss mingw-w64-x86_64-windows-default-manifest
MINGW64_WINDOWS_DEFAULT_MANIFEST_VERSION = 6.4-4
