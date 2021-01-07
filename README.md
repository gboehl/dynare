<a name="logo"/>
<div align="center">
<a href="https://www.dynare.org/" target="_blank">
<img src="https://www.dynare.org/assets/images/logo/dlogo.svg" alt="Dynare Logo"></img>
</a>
</div>

# Dynare

Described on the homepage: <https://www.dynare.org/>

Most users should use the precompiled package available for their OS, also
available via the Dynare homepage: <https://www.dynare.org/download/>.

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
- using the stable source archive of the latest Dynare version from <https://www.dynare.org/download/>
- using a source snapshot of the unstable version, also from <https://www.dynare.org/download/>

Note that if you obtain the source code via git, you will need to install more tools (see below).

The first section of this page gives general instructions, which apply to all platforms. Then some specific platforms are discussed.

**Note:** Here, when we refer to 32-bit or 64-bit, we refer to the type of
MATLAB or Octave installation, not the type of operating system installation.
For example, it is perfectly possible to run a 32-bit MATLAB on a 64-bit
Windows: in that case, instructions for Windows 32-bit should be followed. To
determine the type of your MATLAB/Octave installation, type:
```matlab
>> computer
```
at the MATLAB/Octave prompt. Under MATLAB, if it returns `PCWIN64`, `GLNX64` or
`MACI64`, then it is a 64-bit MATLAB; if it returns `PCWIN`, `MACI` or `GLNX`,
then it is a 32-bit MATLAB. Under Octave, if it returns a string that begins
with `x86_64`, it is a 64-bit Octave; if the strings begins with `i686`, it is
a 32-bit Octave.

**Contents**

1. [**General Instructions**](#general-instructions)
1. [**Debian or Ubuntu**](#debian-or-ubuntu)
1. [**Fedora, CentOS or RHEL**](#fedora-centos-or-rhel)
1. [**Windows**](#windows)
1. [**macOS**](#macos)

## General Instructions

### Prerequisites

A number of tools and libraries are needed in order to recompile everything. You don't necessarily need to install everything, depending on what you want to compile.

- A POSIX compliant shell and an implementation of Make (mandatory)
- The [GNU Compiler Collection](http://gcc.gnu.org/), version 8 or later, with
  gcc, g++ and gfortran (mandatory)
- [MATLAB](https://mathworks.com) (if you want to compile the MEX for MATLAB)
- [GNU Octave](http://www.octave.org) with
  - the development headers (if you want to compile the MEX for Octave)
  - the development libraries corresponding to the [UMFPACK](https://people.engr.tamu.edu/davis/suitesparse.html) packaged with Octave
  - Optionally, the [Control](https://wiki.octave.org/Control_package), [IO](https://wiki.octave.org/IO_package), [Optimization](https://wiki.octave.org/Optimization_package) and [Statistics](https://wiki.octave.org/Statistics_package) package either installed via your package manager or through [Octave Forge](https://wiki.octave.org/Category:Octave_Forge).

- [Boost libraries](http://www.boost.org), version 1.36 or later
- [Bison](http://www.gnu.org/software/bison/), version 3.2 or later (only if you get the source through Git)
- [Flex](http://flex.sourceforge.net/), version 2.5.4 or later (only if you get the source through Git)
- [Autoconf](http://www.gnu.org/software/autoconf/), version 2.62 or later (only if you get the source through Git)
- [Automake](http://www.gnu.org/software/automake/), version 1.11.2 or later (only if you get the source through Git)
- An implementation of BLAS and LAPACK: either [ATLAS](http://math-atlas.sourceforge.net/), [OpenBLAS](http://xianyi.github.com/OpenBLAS/), Netlib ([BLAS](http://www.netlib.org/blas/), [LAPACK](http://www.netlib.org/lapack/)) or [MKL](http://software.intel.com/en-us/intel-mkl/) (only if you want to build Dynare++)
- [MAT File I/O library](http://sourceforge.net/projects/matio/), version 1.5 or later (if you want to compile Markov-Switching code, the estimation DLL, k-order DLL and Dynare++)
- [SLICOT](http://www.slicot.org) (if you want to compile the Kalman steady state DLL)
- [GSL library](http://www.gnu.org/software/gsl/) (if you want to compile Markov-Switching code)
- A decent LaTeX distribution (if you want to compile PDF documentation),
  ideally with Beamer
- For building the reference manual:
  - [Sphinx](http://www.sphinx-doc.org/)
  - [MathJax](https://www.mathjax.org/)
- [Doxygen](http://www.stack.nl/%7Edimitri/doxygen/) (if you want to build Dynare preprocessor source documentation)
- [X-13ARIMA-SEATS Seasonal Adjustment Program](https://www.census.gov/srd/www/x13as/)

### Preparing the sources

If you have downloaded the sources from an official source archive or the source snapshot, just unpack it.

If you want to use Git, do the following from a terminal:
```sh
git clone --recurse-submodules https://git.dynare.org/Dynare/dynare.git
cd dynare
autoreconf -si
```
The last line runs Autoconf and Automake in order to prepare the build environment (this is not necessary if you got the sources from an official source archive or the source snapshot). If you want a certain version (e.g. 4.6) , then add `--single-branch --branch 4.6` to the git clone command.

### Configuring the build tree

Simply launch the configure script from a terminal:
```sh
./configure --with-matlab=<…>
```
where the path to MATLAB is specified.

Some important options:

- `--disable-matlab`: skip the compilation of MEX files for MATLAB
- `--disable-octave`: skip the compilation of MEX files for Octave
- `--disable-doc`: skip the compilation of the documentation (PDF and HTML)

You may need to specify additional options to the configure script, see the output of the `--help` option, and also the platform specific instructions below. If the configuration goes well, the script will tell you which components are correctly configured and will be built. 

Note that it is possible that some MEX files cannot be compiled, due to missing
build dependencies. If you find no way of installing the missing dependencies,
a workaround can be to give up on compiling these MEX files and rather use
slower implementations (in the MATLAB/Octave language) that are available under
the `matlab/missing/mex/` subdirectories. For example, if you fail to compile
the gensylv MEX, you can type the following at the MATLAB/Octave prompt before
running Dynare:
```matlab
addpath <DYNARE_ROOT>/matlab/missing/mex/gensylv
```
(where you need to replace `<DYNARE_ROOT>` with the full path to your Dynare copy).

### Building

Binaries are built with:
```sh
make
```
PDF and HTML documentation are respectively built with:
```sh
make pdf
make html
```

### Check

The Git source comes with unit tests (in the MATLAB functions) and integration tests (under the `tests` subfolder). All the tests can be run with:
```sh
make check
```
in the `tests` subfolder. If Dynare has been compiled against MATLAB and Octave, the tests will be run with both MATLAB and Octave. Depending on the performance of your machine, this can take several hours. It is possible to run the tests only with MATLAB:
```sh
make check-matlab
```
or only with Octave:
```sh
make check-octave
```
Note that running the testsuite with Octave requires the additional packages `pstoedit`, `epstool`, `xfig`, and `gnuplot`. 

A summary of the results is available in `tests/run_test_matlab_output.txt` or `tests/run_test_octave_output.txt`. Often, it does not make sense to run the complete testsuite. For instance, if you modify codes only related to the perfect foresight model solver, you can decide to run only a subset of the integration tests, with:
```sh
make deterministic_simulations
```
This will run all the integration tests in `tests/deterministic_simulations` with MATLAB and Octave. Again, it is possible to do this only with MATLAB:
```sh
make m/deterministic_simulations
```
or with Octave:
```sh
make o/deterministic_simulations
```
Finally if you want to run a single integration test, e.g. `deterministic_simulations/lbj/rbc.mod` with MATLAB:
```sh
make deterministic_simulations/lbj/rbc.m.trs
```
or with Octave:
```sh
make deterministic_simulations/lbj/rbc.o.trs
```
The result of the test (`PASSED` or `FAILED`) will be printed in the terminal, the produced log can be displayed with:
```sh
make deterministic_simulations/lbj/rbc.m.drs
```
or
```sh
make deterministic_simulations/lbj/rbc.o.drs
```
Note that only tests will be executed where the `m.trs/o.trs` does not yet exist. You can run
```sh
make clean
```
in the `tests` folder to delete files that were created by the run of the testsuite. You can also manually delete the desired `m.trs/o.trs` file(s).

## Debian or Ubuntu

All the prerequisites are packaged:

- `build-essential` (for gcc, g++ and make)
- `gfortran`
- `liboctave-dev`
- `libboost-graph-dev`
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
- `texlive-plain-generic`
- `lmodern` (for macroprocessor PDF)
- `python3-sphinx`
- `latexmk`
- `libjs-mathjax`
- `doxygen`
- `x13as`

You can install them all at once with:
```sh
apt install build-essential gfortran liboctave-dev libboost-graph-dev libgsl-dev libmatio-dev libslicot-dev libslicot-pic libsuitesparse-dev flex bison autoconf automake texlive texlive-publishers texlive-latex-extra texlive-fonts-extra texlive-latex-recommended texlive-science texlive-plain-generic lmodern python3-sphinx latexmk libjs-mathjax doxygen x13as
```
If you use MATLAB, we strongly advise to also `apt install matlab-support` and confirm to rename the GCC libraries shipped with MATLAB to avoid possible conflicts with GCC libraries shipped by your distribution.

Tested on
- Debian 10
- Ubuntu 20.04
- Ubuntu 20.10

## Fedora, CentOS or RHEL

Almost all prerequisites are packaged:

- `gcc`, `gcc-c++`, `make`
- `gcc-gfortran`
- `lapack` and `lapack-devel`
- `openblas` and `openblas-devel`
- `boost-devel`
- `gsl-devel`
- `matio-devel`
- `suitesparse-devel`
- `flex`
- `bison`
- `autoconf`
- `automake`
- `redhat-rpm-config`
-  `octave`, `octave-devel`, `octave-statistics`, `octave-io`, `octave-optim`, `octave-control`
- `texlive-scheme-minimal`, `texlive-collection-publishers`, `texlive-collection-latexextra`, `texlive-collection-fontsextra`, `texlive-collection-latexrecommended`, `texlive-collection-science`, `texlive-collection-plaingeneric`, `texlive-lm`
- `python3-sphinx`
- `latexmk`
- `mathjax`
- `doxygen`

You can install them all at once with:
```sh
# Minimal packages (use --disable-doc and --disable-octave flags)
dnf install -y gcc gcc-c++ make gcc-gfortran lapack lapack-devel openblas openblas-devel boost-devel gsl-devel matio-devel suitesparse-devel flex bison autoconf automake redhat-rpm-config
# Octave packages (use --disable-doc flag)
dnf install octave octave-devel octave-statistics octave-io octave-optim octave-control
# Documentation packages
dnf install texlive-scheme-minimal texlive-collection-publishers texlive-collection-latexextra texlive-collection-fontsextra texlive-collection-latexrecommended texlive-collection-science texlive-collection-plaingeneric texlive-lm python3-sphinx latexmk mathjax doxygen
```
In Fedora these are available from the default repositories; whereas for CentOS and RHEL you need to enable the [Extra Packages for Enterprise Linux (EPEL)](https://fedoraproject.org/wiki/EPEL) repository and either the PowerTools repository for CentOS or the CodeReady Linux Builder repository for RHEL:
```sh
yum install https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm
# CentOS 8
dnf config-manager --set-enabled PowerTools
# RHEL 8
ARCH=$( /bin/arch )
subscription-manager repos --enable "codeready-builder-for-rhel-8-${ARCH}-rpms"
```
The documentation packages have slightly different names in CentOS and RHEL, you can also choose to pass the `--disable-doc` flag to your configure script to skip these dependencies.

`Slicot` and `x13as` need to be compiled from source:

```sh
# compile slicot from source and put it into /home/$USER/dynare/slicot/lib/
mkdir -p /home/$USER/dynare/slicot
cd /home/$USER/dynare/slicot
wget https://deb.debian.org/debian/pool/main/s/slicot/slicot_5.0+20101122.orig.tar.gz
tar xf slicot_5.0+20101122.orig.tar.gz
cd slicot-5.0+20101122
make FORTRAN=gfortran OPTS="-O2 -fPIC -fdefault-integer-8" LOADER=gfortran lib
mkdir -p /home/$USER/dynare/slicot/lib
cp slicot.a /home/$USER/dynare/slicot/lib/libslicot64_pic.a #for matlab
cp slicot.a /home/$USER/dynare/slicot/lib/libslicot_pic.a #for octave

# compile x13as from source and put it into /usr/bin/
mkdir -p /home/$USER/dynare/x13as
cd /home/$USER/dynare/x13as
wget https://www.census.gov/ts/x13as/unix/x13assrc_V1.1_B39.tar.gz
tar xf x13assrc_V1.1_B39.tar.gz
sed -i "s|-static| |" makefile.gf # this removes '-static' in the makefile.gf
make -f makefile.gf FFLAGS="-O2 -std=legacy" PROGRAM=x13as
sudo cp x13as /usr/bin/
```

If you use MATLAB, we strongly advise to also rename or exclude the GCC libraries shipped with MATLAB to avoid possible conflicts with GCC libraries shipped by Fedora, see e.g. [Matlab on Fedora 33](https://mutschler.eu/linux/install-guides/fedora-post-install/#matlab) or [MATLAB-ArchWiki](https://wiki.archlinux.org/index.php/MATLAB) for instructions.

Keep in mind to use the `--with-slicot` option to the configure command, e.g.:
```sh
cd /home/$USER/dynare
git clone --recurse-submodules https://git.dynare.org/dynare/dynare.git unstable
cd unstable
autoreconf -si
./configure --with-slicot=/home/$USER/dynare/slicot --with-matlab=/usr/local/MATLAB/R2020b
make -j$(($(nproc)+1)) #rule of thumb: one more than CPUs as shown by e.g. lscpu
```

If your distribution ships an older version of `bison`, compile it from source and append it *temporarily* to your path before calling the configure script:
```sh
bison --version # bison (GNU Bison) 3.0.4
mkdir -p /home/$USER/dynare/bison
cd /home/$USER/dynare/bison
wget http://ftp.gnu.org/gnu/bison/bison-3.6.4.tar.gz #change the version number accordingly
tar xf bison-3.6.4.tar.gz
cd bison-3.6.4
./configure --prefix=/home/$USER/dynare/bison
make
make install
export PATH=/home/$USER/dynare/bison/bin:$PATH
bison --version # bison (GNU Bison) 3.6.4
```
Now configure dynare as above.

Tested on
- CentOS 8
- Fedora Workstation 32
- Fedora Workstation 33
- Red Hat Enterprise Linux 8

## Windows

- Install [MSYS2](http://www.msys2.org) (pick the 64-bit version, unless you
  have a 32-bit Windows, in which case see below)
- Run a MSYS MinGW 64-bit shell
- Update the system:
```sh
pacman -Syu
```
  You may be asked to close the window at the end of the
  first upgrade batch, in which case you should rerun the upgrade in a new
  window to complete the upgrade.
- Install all needed dependencies:
```sh
pacman -S git autoconf automake-wrapper bison flex make tar texinfo mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-boost mingw-w64-x86_64-gsl mingw-w64-x86_64-matio mingw-w64-x86_64-openblas
```
- Compile and install SLICOT, needed for the `kalman_steady_state` MEX file
```sh
wget https://deb.debian.org/debian/pool/main/s/slicot/slicot_5.0+20101122.orig.tar.gz
tar xf slicot_5.0+20101122.orig.tar.gz
cd slicot-5.0+20101122
make FORTRAN=gfortran OPTS="-O2 -fno-underscoring -fdefault-integer-8" LOADER=gfortran lib
mkdir -p /usr/local/lib
cp slicot.a /usr/local/lib/libslicot64_pic.a
cd ..
```
- Prepare the Dynare sources, either by unpacking the source tarball, or with:
```sh
git clone --recurse-submodules https://git.dynare.org/Dynare/dynare.git
cd dynare
autoreconf -si
```
- Configure Dynare from the source directory:
```sh
./configure --with-slicot=/usr/local --with-matlab=<…> --disable-octave --disable-doc
```
where the path of MATLAB is specified. Note that you should use
the MSYS2 notation and not put spaces in the MATLAB path, so you probably want
to use something like `/c/Progra~1/MATLAB/…`. Alternatively, if your filesystem
does not have short filenames (8dot3), then you can run `mkdir -p
/usr/local/MATLAB && mount c:/Program\ Files/MATLAB /usr/local/MATLAB`, and
then pass `/usr/local/MATLAB/…` as MATLAB path to the configure script.
- Compile:
```sh
make
```
- Run the testsuite:
```sh
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

Preparatory work:

- Install the Xcode Command Line Tools. Open Terminal.app and type:
```sh
xcode-select --install
```
- Install [Homebrew](https://brew.sh/) by following the instructions on their website
- Install [MacTeX](http://www.tug.org/mactex/index.html). Alternatively, if you
  don’t want to install MacTeX, you should pass the `--disable-doc` flag to the
  `configure` command below.
- Install required Homebrew packages. Open Terminal.app and type:
```sh
brew install automake bison flex boost gcc gsl libmatio veclibfort octave sphinx-doc wget
brew link --force sphinx-doc
```
- Compile and install SLICOT, needed for the `kalman_steady_state` MEX file.
Still from Terminal.app:
```sh
wget https://deb.debian.org/debian/pool/main/s/slicot/slicot_5.0+20101122.orig.tar.gz
tar xf slicot_5.0+20101122.orig.tar.gz
cd slicot-5.0+20101122
make -j$(nproc) FORTRAN=gfortran OPTS="-O2" LOADER=gfortran lib
cp slicot.a /usr/local/lib/libslicot_pic.a
make clean
make -j$(nproc) FORTRAN=gfortran OPTS="-O2 -fdefault-integer-8" LOADER=gfortran lib
cp slicot.a /usr/local/lib/libslicot64_pic.a
cd ..
```

The following commands will download the Dynare source code and compile
it. They should be entered at the command prompt in Terminal.app from the
folder where you want Dynare installed.

- Prepare the Dynare sources:
```sh
git clone --recurse-submodules https://git.dynare.org/Dynare/dynare.git`
cd dynare
autoreconf -si
```
- Configure Dynare from the source directory:
```sh
./configure --with-matlab=<…> CC=gcc-10 CXX=g++-10 CPPFLAGS=-I/usr/local/include LDFLAGS=-L/usr/local/lib` LEX=/usr/local/opt/flex/bin/flex YACC=/usr/local/opt/bison/bin/bison'
```
where the path to MATLAB is specified, typically of the form
`/Applications/MATLAB_R2020b.app`. If you don’t have MATLAB, simply replace
`--with-matlab=<…>` in the above command by `--disable-matlab`.
- Compile:
```sh
make -j$(nproc)
```
