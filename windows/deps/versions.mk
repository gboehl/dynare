SLICOT_VERSION = 5.0+20101122
X13AS_VERSION = 1.1_B39

OCTAVE_VERSION = 6.4.0
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
MINGW64_BOOST_VERSION = 1.78.0-3

# pacman -Ss mingw-w64-x86_64-gsl
MINGW64_GSL_VERSION = 2.7.1-1

# pacman -Ss mingw-w64-x86_64-openblas
MINGW64_OPENBLAS_VERSION = 0.3.19-1

# pacman -Ss mingw-w64-x86_64-matio
MINGW64_MATIO_VERSION = 1.5.21-2

# Dependency of matio (and of the MinGW compiler)
# pacman -Ss mingw-w64-x86_64-zlib
MINGW64_ZLIB_VERSION = 1.2.11-9

# Dependency of matio
# pacman -Ss mingw-w64-x86_64-hdf5
MINGW64_HDF5_VERSION = 1.12.1-2

# Dependency of HDF5 (provides szip library)
# pacman -Ss mingw-w64-x86_64-libaec
MINGW64_LIBAEC_VERSION = 1.0.6-2

## MinGW packages for the embedded compiler

# pacman -Ss mingw-w64-x86_64-gcc$
MINGW64_GCC_VERSION = 11.2.0-9

# pacman -Ss mingw-w64-x86_64-gmp
MINGW64_GMP_VERSION = 6.2.1-3

# pacman -Ss mingw-w64-x86_64-binutils
MINGW64_BINUTILS_VERSION = 2.38-1

# pacman -Ss mingw-w64-x86_64-headers-git
MINGW64_HEADERS_VERSION = 9.0.0.6373.5be8fcd83-2

# pacman -Ss mingw-w64-x86_64-crt-git
MINGW64_CRT_VERSION = 9.0.0.6373.5be8fcd83-2

# pacman -Ss mingw-w64-x86_64-winpthreads-git
MINGW64_WINPTHREADS_VERSION = 9.0.0.6373.5be8fcd83-2

# pacman -Ss mingw-w64-x86_64-zstd
MINGW64_ZSTD_VERSION = 1.5.2-1
