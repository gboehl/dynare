<a name="logo"/>
<div align="center">
<a href="http://www.dynare.org/" target="_blank">
<img src="http://www.dynare.org/img/dynare.png" alt="Dynare Logo"></img>
</a>
</div>

# Dynare

Described on the homepage: <http://www.dynare.org/>

Most users should use the precompiled package available for your OS, also
available via the Dynare homepage: <http://www.dynare.org/download/dynare-stable>.

# Contributions

To contribute to Dynare and participate in the Dynare community, please see: [CONTRIBUTING.md](https://git.dynare.org/Dynare/dynare/blob/master/CONTRIBUTING.md)

# License

Most of the source files are covered by the GNU General Public Licence version
3 or later (there are some exceptions to this, see [license.txt](license.txt) in
Dynare distribution for specifics).

# Building Dynare From Source

Here, we explain how to build from source:
- Dynare, including preprocessor and MEX files for MATLAB and Octave
- Dynare++
- all the associated documentation (PDF and HTML)

This source can be retrieved in three forms:
- via git, at <https://git.dynare.org/Dynare/dynare.git>
- using the stable source archive of the latest Dynare version (currently 4.4) from <http://www.dynare.org/download/dynare-stable/>
- using a source snapshot of the unstable version, from <http://www.dynare.org/download/dynare-unstable/source-snapshot>

Note that if you obtain the source code via git, you will need to install more tools (see below).

The first section of this page gives general instructions, which apply to all platforms. Then some specific platforms are discussed.

**Note:** Here, when we refer to 32-bit or 64-bit, we refer to the type of MATLAB installation, not the type of Windows installation. It is perfectly possible to run a 32-bit MATLAB on a 64-bit Windows: in that case, instructions for Windows 32-bit should be followed. To determine the type of your MATLAB installation, type:
```matlab
>> computer
```
at the MATLAB prompt: if it returns `PCWIN64`, `GLNX64` or `MACI64`, then you
have a 64-bit MATLAB; if it returns `PCWIN`, `MACI` or `GLNX`, then you have a
32-bit MATLAB.

**Contents**

1. [**General Instructions**](#general-instructions)
1. [**Debian or Ubuntu**](#debian-or-ubuntu)
1. [**Windows**](#windows)
1. [**macOS**](#macos)

## General Instructions

### Prerequisites

A number of tools and libraries are needed in order to recompile everything. You don't necessarily need to install everything, depending on what you want to compile.

- A POSIX compliant shell and an implementation of Make (mandatory)
- The [GNU Compiler Collection](http://gcc.gnu.org/), version 6 or later, with
  gcc, g++ and gfortran (mandatory)
- MATLAB (if you want to compile the MEX for MATLAB)
- [GNU Octave](http://www.octave.org), with the development headers (if you
  want to compile the MEX for Octave)
- [Boost libraries](http://www.boost.org), version 1.36 or later (with the filesystem library compiled)
- [Bison](http://www.gnu.org/software/bison/), version 3.0 or later (only if you get the source through Git)
- [Flex](http://flex.sourceforge.net/), version 2.5.4 or later (only if you get the source through Git)
- [Autoconf](http://www.gnu.org/software/autoconf/), version 2.62 or later (only if you get the source through Git) (see [Installing an updated version of Autoconf in your own directory, in GNU/Linux](http://www.dynare.org/DynareWiki/AutoMake))
- [Automake](http://www.gnu.org/software/automake/), version 1.11.2 or later (only if you get the source through Git) (see [Installing an updated version of AutoMake in your own directory, in GNU/Linux](http://www.dynare.org/DynareWiki/AutoMake))
- An implementation of BLAS and LAPACK: either [ATLAS](http://math-atlas.sourceforge.net/), [OpenBLAS](http://xianyi.github.com/OpenBLAS/), Netlib ([BLAS](http://www.netlib.org/blas/), [LAPACK](http://www.netlib.org/lapack/)) or [MKL](http://software.intel.com/en-us/intel-mkl/) (only if you want to build Dynare++)
- [MAT File I/O library](http://sourceforge.net/projects/matio/), version 1.5 or later (if you want to compile Markov-Switching code, the estimation DLL, k-order DLL and Dynare++)
- [SLICOT](http://www.slicot.org) (if you want to compile the Kalman steady state DLL)
- [GSL library](http://www.gnu.org/software/gsl/) (if you want to compile Markov-Switching code)
- A decent LaTeX distribution (if you want to compile PDF documentation),
  ideally with Beamer
- For building the reference manual:
  - [GNU Texinfo](http://www.gnu.org/software/texinfo/)
  - [Latex2HTML](http://www.latex2html.org), if you want nice mathematical formulas in HTML output
  - [Doxygen](http://www.stack.nl/%7Edimitri/doxygen/) (if you want to build Dynare preprocessor source documentation)
- For Octave, the development libraries corresponding to the UMFPACK packaged with Octave

### Preparing the sources

If you have downloaded the sources from an official source archive or the source snapshot, just unpack it.

If you want to use Git, do the following from a terminal:

    git clone --recursive https://git.dynare.org/Dynare/dynare.git
    cd dynare
    autoreconf -si

The last line runs Autoconf and Automake in order to prepare the build environment (this is not necessary if you got the sources from an official source archive or the source snapshot).

### Configuring the build tree

Simply launch the configure script from a terminal:
```
./configure
```
If you have MATLAB, you need to indicate both the MATLAB location and version. For example, on GNU/Linux:
```
./configure --with-matlab=/usr/local/MATLAB/R2013a MATLAB_VERSION=8.1
```
Note that the MATLAB version can also be specified via the MATLAB family product release (R2009a, R2008b, ...).

Alternatively, you can disable the compilation of MEX files for MATLAB with the `--disable-matlab` flag, and MEX files for Octave with `--disable-octave`.

You may need to specify additional options to the configure script, see the platform specific instructions below.

Note that if you don't want to compile the C/C++ programs with debugging information, you can specify the `CFLAGS` and `CXXFLAGS` variables to the configure script, such as:
```
./configure CFLAGS="-O3" CXXFLAGS="-O3"
```
To remove debugging information for MATLAB MEX functions, the analagous call would be:
```
./configure MATLAB_MEX_CFLAGS="-O3" MATLAB_MEX_CXXFLAGS="-O3"
```

If you want to give a try to the parallelized versions of some mex files (`A_times_B_kronecker_C` and `sparse_hessian_times_B_kronecker_C` used to get the reduced form of the second order approximation of the model) you can add the `--enable-openmp` flag, for instance:
```
./configure --with-matlab=/usr/local/MATLAB/R2013a MATLAB_VERSION=8.1 --enable-openmp
```
If the configuration goes well, the script will tell you which components are correctly configured and will be built.

### Building

Binaries and Info documentation are built with:
```
make
```
PDF and HTML documentation are respectively built with:
```
make pdf
make html
```
The testsuites can be run with:
```
make check
```

Note that running the testsuite with Octave requires the additional packages
`pstoedit`, `epstool`, `xfig`, and `gnuplot`.

### Check

The Git source comes with unit tests (in the MATLAB functions) and integration tests (under the `tests` subfolder). All the tests can be run with:
```
make check
```
In the `tests` subfolder. If Dynare has been compiled against MATLAB and Octave, the tests will be run with MATLAB and Octave. Depending on
your PC, this can take several hours. It is possible to run the tests only with MATLAB:
```
make check-matlab
```
or only with Octave:
```
make check-octave
```
A summary of the results is available in `tests/run_test_matlab_output.txt` or `tests/run_test_octave_output.txt`. Often, it does not make sense
to run the complete testsuite. For instance, if you modify codes only related to the perfect foresight model solver, you can decide to run only a
subset of the integration tests, with:
```
make deterministic_simulations
```
This will run all the integration tests in `tests/deterministic_simulations` with MATLAB and Octave. Again, it is possible to do this only with MATLAB:
```
make m/deterministic_simulations
```
or with Octave:
```
make o/deterministic_simulations
```
Finally if you want to run a single integration test, e.g. `deterministic_simulations/lbj/rbc.mod` with MATLAB:
```
make deterministic_simulations/lbj/rbc.m.trs
```
or with Octave:
```
make deterministic_simulations/lbj/rbc.o.trs
```
The result of the test (`PASSED` or `FAILED`) will be printed in the terminal, the produced log can be displayed with:
```
make deterministic_simulations/lbj/rbc.m.drs
```
or
```
make deterministic_simulations/lbj/rbc.o.drs
```
Note that only tests will be executed where the `m.trs/o.trs` does not yet exist. You can run
```
make clean
```
in the `tests` folder to delete files that were created by the run of the testsuite. You can also manually delete the desired `m.trs/o.trs` file(s).

## Debian or Ubuntu

All the prerequisites are packaged:

- `build-essential` (for gcc, g++ and make)
- `gfortran`
- `liboctave-dev`
- `libboost-graph-dev` and `libboost-filesystem-dev`
- `libgsl-dev`
- `libmatio-dev`
- `libslicot-dev` and `libslicot-pic`
- `libsuitesparse-dev`
- `flex`
- `bison`
- `autoconf`
- `automake`
- `texlive`
- `texlive-publishers` (for Econometrica bibliographic style)
- `texlive-latex-extra` (for fullpage.sty)
- `texlive-fonts-extra` (for ccicons)
- `texlive-latex-recommended`
- `texlive-science` (for amstex)
- `texinfo`
- `lmodern` (for macroprocessor PDF)
- `latex2html`
- `doxygen`

You can install them all at once with:
```
apt install build-essential gfortran liboctave-dev libboost-graph-dev libboost-filesystem-dev libgsl-dev libmatio-dev libslicot-dev libslicot-pic libsuitesparse-dev flex bison autoconf automake texlive texlive-publishers texlive-extra-utils texlive-formats-extra texlive-latex-extra texlive-fonts-extra texlive-latex-recommended texlive-science texinfo lmodern latex2html doxygen
```

## Windows

- Install [MSYS2](http://www.msys2.org) (pick the 64-bit version)
- Run a MSYS MinGW 64-bit shell
- Update the system:
```
pacman -Syu
```
  You may be asked to close the window at the end of the
  first upgrade batch, in which case you should rerun the upgrade in a new
  window to complete the upgrade.
- Install all needed dependencies:
```
pacman -S git autoconf automake-wrapper bison flex make tar texinfo mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-boost mingw-w64-x86_64-gsl mingw-w64-x86_64-matio mingw-w64-x86_64-openblas
```
- **(Optional)** compile and install SLICOT, needed for the `kalman_steady_state`
  MEX file
```
wget https://deb.debian.org/debian/pool/main/s/slicot/slicot_5.0+20101122.orig.tar.gz
tar xf slicot_5.0+20101122.orig.tar.gz
cd slicot-5.0+20101122
make FORTRAN=gfortran OPTS="-O2 -fno-underscoring -fdefault-integer-8" LOADER=gfortran slicot.a
mkdir -p /usr/local/lib
cp slicot.a /usr/local/lib/libslicot64_pic.a
cd ..
```
- Clone and prepare the Dynare sources:
```
git clone --recursive https://git.dynare.org/Dynare/dynare.git
cd dynare
autoreconf -si
```
- Configure Dynare:
```
./configure --with-boost-system=boost_system-mt --with-boost-filesystem=boost_filesystem-mt --with-slicot=/usr/local --with-matlab=<…> MATLAB_VERSION=<…> --disable-octave
```
where the path and version of MATLAB are specified. Note that you should use
the MSYS2 notation and not put spaces in the MATLAB path, so you probably want
to use something like `/c/Progra~1/MATLAB/…`.
- Compile:
```
make
```
- Run the testsuite:
```
make -C tests check-matlab
```

**Note:** The above assumes that you have a 64-bit version of MATLAB. It can be
adapted to a 32-bit MATLAB with the following modifications:

- run the MSYS MinGW 32-bit shell
- replace `x86_64` by `i686` in packages names on the `pacman` command-line
- for SLICOT, remove the `-fdefault-integer-8` option, and instead copy the
  library into `/usr/local/lib/libslicot_pic.a`

**Note:** Compiling the MEX files for Octave and the documentation under MSYS2 is
currently not supported.

## macOS

To simply use a snapshot of Dynare, you have two choices. On MATLAB, you can
use the [snapshot build](http://www.dynare.org/snapshot/macosx/) provided by
Dynare. On Octave, you can simply install [Homebrew](https://brew.sh/) and run
```brew install dynare --HEAD``` (See the Install Dynare (unstable) section of
[this webpage](http://www.dynare.org/DynareWiki/InstallOnMacOSX) for more
details).

If you do not wish to use the snapshots provided by Dynare or Homebrew, follow
the directions below to build Dynare on your local machine.

Preparatory work:

- Install the Xcode Command Line Tools:
    - Open Terminal.app and type `xcode-select --install`
- Install [Homebrew](https://brew.sh/) by following the instructions on their website

The following commands will install the programs that Dynare needs to
compile. They should be entered at the command prompt in Terminal.app.

- `brew install automake bison flex boost fftw gcc gsl hdf5 libmatio metis veclibfort`
- **(Optional)** To compile Dynare mex files for use on Octave:
    - `brew install octave`
    - `brew install suite-sparse`
- **(Optional)** To compile Dynare documentation
     - Install the latest version of [MacTeX](http://www.tug.org/mactex/), deselecting the option to install Ghostscript
     - `brew install doxygen latex2html`

The following commands will download the Dynare source code and compile
it. They should be entered at the command prompt in Terminal.app from the
folder where you want Dynare installed.

- `git clone https://git.dynare.org/Dynare/dynare.git`
- `cd dynare`
- `PATH="/usr/local/opt/bison/bin:/usr/local/opt/flex/bin:$PATH"`
- `autoreconf -si`
- `./configure --disable-octave --with-matlab=/Applications/MATLAB_R2017b.app MATLAB_VERSION=R2017b`, adjusting the MATLAB path and version to accord with your local installation. If you don't have MATLAB, simply type `./configure --disable-octave`
- `make -j`
- **(Optional)** To then build mex files for Octave, run
     - `cd mex/build/octave`
     - `./configure CXXFLAGS="-std=c++0x"`
     - `make -j`
