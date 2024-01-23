SLICOT_VERSION = 5.8+20230223.git401037e
X13AS_VERSION = 1-1-b60

OCTAVE_VERSION = 8.4.0
OCTAVE_W64_BUILD =

MATLAB64_VERSION = 20231122


### MSYS2 packages
# Determine the versions by:
# - first running: pacman -Sy
# - and then with appropriate queries using: pacman -Ss <regex>
# Dependencies can be determined using: pacman -Si <pkg>
# File lists can be determined using: pacman -Fl <pkg>
# The same information can be gathered from: https://packages.msys2.org/search

## Build dependencies

# pacman -Ss mingw-w64-x86_64-boost
MINGW64_BOOST_VERSION = 1.84.0-1

# pacman -Ss mingw-w64-x86_64-gsl
MINGW64_GSL_VERSION = 2.7.1-2

# pacman -Ss mingw-w64-x86_64-matio
MINGW64_MATIO_VERSION = 1.5.24-1

# Dependency of matio (and of the MinGW compiler)
# pacman -Ss mingw-w64-x86_64-zlib
MINGW64_ZLIB_VERSION = 1.3-1

# Dependency of matio
# pacman -Ss mingw-w64-x86_64-hdf5
MINGW64_HDF5_VERSION = 1.14.2-3

# Dependency of HDF5 (provides szip library)
# pacman -Ss mingw-w64-x86_64-libaec
MINGW64_LIBAEC_VERSION = 1.0.6-2

# Dependency of HDF5
# pacman -Ss mingw-w64-x86_64-openssl
MINGW64_OPENSSL_VERSION = 3.2.0-1

# Dependency of HDF5
# pacman -Ss mingw-w64-x86_64-curl
MINGW64_CURL_VERSION = 8.5.0-1

# Dependency of curl (and of the MinGW compiler)
# pacman -Ss mingw-w64-x86_64-zstd
MINGW64_ZSTD_VERSION = 1.5.5-1

# Dependency of curl
# pacman -Ss mingw-w64-x86_64-brotli
MINGW64_BROTLI_VERSION = 1.1.0-1

# Dependency of curl
# pacman -Ss mingw-w64-x86_64-libpsl
MINGW64_LIBPSL_VERSION = 0.21.2-4

# Dependency of curl and of libpsl
# pacman -Ss mingw-w64-x86_64-libidn2
MINGW64_LIBIDN2_VERSION = 2.3.4-1

# Dependency of curl
# pacman -Ss mingw-w64-x86_64-libssh2
MINGW64_LIBSSH2_VERSION = 1.11.0-2

# Dependency of curl
# pacman -Ss mingw-w64-x86_64-nghttp2
MINGW64_NGHTTP2_VERSION = 1.58.0-1

# Dependency of libpsl and libunistring (and of the MinGW compiler)
# pacman -Ss mingw-w64-x86_64-libiconv
MINGW64_LIBICONV_VERSION = 1.17-3

# Dependency of libpsl and libidn2
# pacman -Ss mingw-w64-x86_64-libunistring
MINGW64_LIBUNISTRING_VERSION = 1.1-1

## MinGW packages for the embedded compiler

# pacman -Ss mingw-w64-x86_64-gcc$
MINGW64_GCC_VERSION = 13.2.0-3

# pacman -Ss mingw-w64-x86_64-gmp
MINGW64_GMP_VERSION = 6.3.0-2

# pacman -Ss mingw-w64-x86_64-binutils
MINGW64_BINUTILS_VERSION = 2.41-3

# pacman -Ss mingw-w64-x86_64-headers-git
MINGW64_HEADERS_VERSION = 11.0.0.r442.ga27e7b27e-1

# pacman -Ss mingw-w64-x86_64-crt-git
MINGW64_CRT_VERSION = 11.0.0.r442.ga27e7b27e-1

# pacman -Ss mingw-w64-x86_64-winpthreads-git
MINGW64_WINPTHREADS_VERSION = 11.0.0.r442.ga27e7b27e-1

# pacman -Ss mingw-w64-x86_64-isl
MINGW64_ISL_VERSION = 0.26-1

# pacman -Ss mingw-w64-x86_64-mpc
MINGW64_MPC_VERSION = 1.3.1-2

# pacman -Ss mingw-w64-x86_64-mpfr
MINGW64_MPFR_VERSION = 4.2.1-2

# pacman -Ss mingw-w64-x86_64-windows-default-manifest
MINGW64_WINDOWS_DEFAULT_MANIFEST_VERSION = 6.4-4
