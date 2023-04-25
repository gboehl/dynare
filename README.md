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
at the MATLAB/Octave prompt. Under MATLAB, if it returns `PCWIN64`, `GLNX64`,
`MACI64` or `MACA64` then it is a 64-bit MATLAB; if it returns `PCWIN`, `MACI` or `GLNX`,
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
- The [GNU Compiler Collection](https://gcc.gnu.org/), version 10 or later, with
  gcc, g++ and gfortran (mandatory)
- [MATLAB](https://mathworks.com) (if you want to compile the MEX for MATLAB)
- [GNU Octave](https://www.octave.org) with
  - the development headers (if you want to compile the MEX for Octave)
  - the development libraries corresponding to the [UMFPACK](https://people.engr.tamu.edu/davis/suitesparse.html) packaged with Octave
  - Optionally, the [Control](https://octave.sourceforge.io/control/), [IO](https://octave.sourceforge.io/io/), [Optimization](https://octave.sourceforge.io/optim/) and [Statistics](https://octave.sourceforge.io/statistics/) package either installed via your package manager or through [Octave Forge](https://octave.sourceforge.io/).
- [Boost libraries](https://www.boost.org), version 1.36 or later
- [Bison](https://www.gnu.org/software/bison/), version 3.2 or later (only if you get the source through Git)
- [Flex](https://github.com/westes/flex), version 2.5.4 or later (only if you get the source through Git)
- [Autoconf](https://www.gnu.org/software/autoconf/), version 2.62 or later (only if you get the source through Git)
- [Automake](https://www.gnu.org/software/automake/), version 1.11.2 or later (only if you get the source through Git)
- [MAT File I/O library](https://sourceforge.net/projects/matio/), version 1.5 or later (if you want to compile Markov-Switching SBVAR code)
- [SLICOT](http://www.slicot.org) (if you want to compile the Kalman steady state DLL)
- [GSL library](https://www.gnu.org/software/gsl/) (if you want to compile Markov-Switching SBVAR code)
- A decent LaTeX distribution (if you want to compile PDF documentation),
  ideally with Beamer
- For building the reference manual:
  - [Sphinx](https://www.sphinx-doc.org/)
  - [MathJax](https://www.mathjax.org/)
- [Doxygen](https://www.doxygen.nl) (if you want to build Dynare preprocessor source documentation)
- [X-13ARIMA-SEATS Seasonal Adjustment Program](https://www.census.gov/data/software/x13as.html)

### Preparing the sources

If you have downloaded the sources from an official source archive or the source snapshot, just unpack it.

If you want to use Git, do the following from a terminal (note that you must
have the [Git LFS](https://git-lfs.github.com/) extension installed):
```sh
git clone --recurse-submodules https://git.dynare.org/Dynare/dynare.git
cd dynare
autoreconf -si
```
The last line runs Autoconf and Automake in order to prepare the build environment (this is not necessary if you got the sources from an official source archive or the source snapshot). If you want a certain version (e.g. 5.x) , then add `--single-branch --branch 5.x` to the git clone command.

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
- `flex` and `libfl-dev`
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
- `tex-gyre`
- `latexmk`
- `libjs-mathjax`
- `doxygen`
- `x13as`

You can install them all at once with:
```sh
apt install build-essential gfortran liboctave-dev libboost-graph-dev libgsl-dev libmatio-dev libslicot-dev libslicot-pic libsuitesparse-dev flex libfl-dev bison autoconf automake texlive texlive-publishers texlive-latex-extra texlive-fonts-extra texlive-latex-recommended texlive-science texlive-plain-generic lmodern python3-sphinx tex-gyre latexmk libjs-mathjax doxygen x13as
```
If you use MATLAB, we strongly advise to also `apt install matlab-support` and confirm to rename the GCC libraries shipped with MATLAB to avoid possible conflicts with GCC libraries shipped by your distribution.

Tested on
- Debian “Buster” 10
- Debian “Bullseye” 11
- Ubuntu 20.04
- Ubuntu 20.10

## Fedora, CentOS or RHEL

Almost all prerequisites are packaged:

- `gcc`, `gcc-c++`, `make`
- `gcc-gfortran`
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
dnf install -y gcc gcc-c++ make gcc-gfortran boost-devel gsl-devel matio-devel suitesparse-devel flex bison autoconf automake redhat-rpm-config
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
wget https://www2.census.gov/software/x-13arima-seats/x13as/unix-linux/program-archives/x13as_asciisrc-v1-1-b59.tar.gz
tar xf x13as_asciisrc-v1-1-b59.tar.gz
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

- Install [MSYS2](http://www.msys2.org)
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
pacman -S git autoconf automake-wrapper bison flex make tar texinfo mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-boost mingw-w64-x86_64-gsl mingw-w64-x86_64-matio
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

Dynare supports both Intel and Apple Silicon chips and is compiled from source
using a [Homebrew](https://brew.sh/) toolchain. If you have an Apple silicon processor
(*M1/M2 PRO/MAX/ULTRA*), you can compile Dynare both for Intel's `x86_64` (using Rosetta 2)
as well as Apple's native `arm64` platform by using the corresponding Homebrew packages.
If you have an Intel chip you can only compile for `x86_64`.

You can check the platform of your current Homebrew installation by e.g. running
`which brew` which should point to `/opt/homebrew/bin/brew` for `arm64` and to
`/usr/local/bin/brew` for `x86_64` systems. In the steps below, we
create a temporary environmental variable `BREWDIR` to ensure that the correct packages are used.

The following commands install all requirements and Dynare from source.
They should be entered at the command prompt in Terminal.app.

### Preparatory work

- Install the Xcode Command Line Tools:
```sh
xcode-select --install
```

- Install Rosetta 2 (Apple Silicon only):
```sh
softwareupdate --install-rosetta --agree-to-license
```

- Install [Homebrew](https://brew.sh/):
Create environmental variables for which platform you want to compile for, i.e. either `arm64` or `x86_64`:

For `arm64` run the following commands:
```sh
export ARCH=arm64
export BREWDIR=/opt/homebrew
```

For `x86_64` run the following commands:
```sh
export ARCH=x86_64
export BREWDIR=/usr/local
```

Install Homebrew using the environmental variables:
```sh
arch -$ARCH /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
The prefix `arch -arm64` or `arch -x86_64` makes sure that you are installing the correct packages.
Don't forget to run the displayed commands (**Next steps**) in the terminal to add Homebrew to your PATH.

If you have both Homebrew installations installed, make sure that you are accessing the correct packages by temporarily (!) prepending it to the path:
```sh
export PATH="$BREWDIR/bin:$PATH"
```


- Install required Homebrew packages and link sphinx-doc:
```sh
arch -$ARCH $BREWDIR/bin/brew install automake bison flex boost gcc gsl libmatio veclibfort octave sphinx-doc docutils wget
```
If you want to compile the documentation, you need to link sphinx-doc:
```sh
arch -$ARCH $BREWDIR/bin/brew link --force sphinx-doc
```
Otherwise add `--disable-doc` flag to the `configure` command below, if you want to skip this.

- Install [MacTeX](http://www.tug.org/mactex/index.html) using the universal installer. MacTeX runs natively on both ARM and Intel machines. On Apple Silicon, it is advised to symlink `pdflatex` and `bibtex` into `/usr/local/bin`:
```sh
sudo ln -s /Library/TeX/texbin/pdflatex /usr/local/bin/pdflatex
sudo ln -s /Library/TeX/texbin/bibtex /usr/local/bin/bibtex
```
Alternatively, if you don’t want to install MacTeX, you should pass the `--disable-doc` flag to the `configure` command below.

- Install MATLAB and additional toolboxes.
We recommend, but don't require, the following: Optimization, Global Optimization, Statistics and Machine Learning, Econometrics, and Control System.
For Apple Silicon: MATLAB offers a native Apple silicon version of R2022b as an open beta. You can sign up and install it (including a suitable Java 8 JRE, e.g. Amazon Corretto 8) using
[the official instructions](https://de.mathworks.com/support/apple-silicon-r2022b-beta.html).
Unfortunately, this version does not support yet the Optimization, Global Optimization and Econometrics toolboxes.
If you need these, please run the Intel version (under Rosetta 2) instead.
Don't forget to run MATLAB at least once to make sure you have a valid license.

- Create a folder for Dynare and its dependencies
```sh
export DYNAREDIR=$HOME/dynare
```

- Compile and install SLICOT, needed for e.g. the `kalman_steady_state` MEX file.
```sh
mkdir -p $DYNAREDIR/slicot/lib
cd $DYNAREDIR/slicot
curl -O https://deb.debian.org/debian/pool/main/s/slicot/slicot_5.0+20101122.orig.tar.gz
tar xf slicot_5.0+20101122.orig.tar.gz
cd slicot-5.0+20101122
make -j$(sysctl -n hw.ncpu) FORTRAN=$BREWDIR/bin/gfortran OPTS="-O2" LOADER=gfortran lib
cp slicot.a $DYNAREDIR/slicot/lib/libslicot_pic.a
make clean
make -j$(sysctl -n hw.ncpu) FORTRAN=$BREWDIR/bin/gfortran OPTS="-O2 -fdefault-integer-8" LOADER=gfortran lib
cp slicot.a $DYNAREDIR/slicot/lib/libslicot64_pic.a
cd $HOME/dynare
```

- Compile and install the X-13ARIMA-SEATS Seasonal Adjustment Program
```sh
mkdir -p $DYNAREDIR/x13as
cd $DYNAREDIR/x13as
curl -O https://www2.census.gov/software/x-13arima-seats/x13as/unix-linux/program-archives/x13as_asciisrc-v1-1-b59.tar.gz
tar xf x13as_asciisrc-v1-1-b59.tar.gz
sed -i '' 's/-static//g' makefile.gf
make -j$(sysctl -n hw.ncpu) -f makefile.gf FC=$BREWDIR/bin/gfortran LINKER=$BREWDIR/bin/gcc-12 FFLAGS="-O2 -std=legacy" LDFLAGS=-static-libgcc LIBS="$BREWDIR/lib/gcc/current/libgfortran.a /$BREWDIR/lib/gcc/current/libquadmath.a" PROGRAM=x13as
mkdir -p $HOME/.local/bin
cp x13as $HOME/.local/bin/x13as
cd ;
x13as
```

### Compile Dynare from source
The following commands will download the Dynare source code and compile
it. They should be entered at the command prompt in Terminal.app from the
folder where you want Dynare installed.

- Prepare the Dynare sources for the unstable version:
```sh
mkdir -p $DYNAREDIR/unstable
git clone --recurse-submodules https://git.dynare.org/Dynare/dynare.git $DYNAREDIR/unstable
cd $DYNAREDIR/unstable
arch -$ARCH autoreconf -si
```
You can also choose a specific version of Dynare by checking out the corresponding branch or a specific tag with git.

- Configure Dynare from the source directory:
```sh
arch -$ARCH ./configure --prefix=$BREWDIR CC=gcc-12 CXX=g++-12 CPPFLAGS=-I$BREWDIR/include LDFLAGS=-L$BREWDIR/lib LEX=$BREWDIR/opt/flex/bin/flex YACC=$BREWDIR/opt/bison/bin/bison --with-slicot=$DYNAREDIR/slicot --with-matlab=/Applications/MATLAB_R2022b_Beta.app
```
where you need to adapt the path to MATLAB. If you don’t have MATLAB, simply replace `--with-matlab=<…>` by `--disable-matlab`.
Similarly, if you don't want to compile for Octave add a `--disable-octave` flag.
Check the output of the command whether Dynare is configured for building everything except the internal docs of Dynare and M2HTML.

- Compile:
```sh
arch -$ARCH make -j$(sysctl -n hw.ncpu)
```
If no errors occured, you are done. Dynare is now ready to use.
If you additionally want to compile the documentation run:
```sh
arch -$ARCH make pdf
arch -$ARCH make html
```

### Optional: pass the full PATH to MATLAB to run system commands
If you start MATLAB from a terminal, you will get the PATH inherited from the shell.
However, when you click on the application icon in macOS, you are not running at the terminal level:
the program is run by launcher, which does not go through a shell login session.
In other words, you get the system default PATH which includes `/usr/bin:/bin:/usr/sbin:/sbin`, but not `/usr/local/bin` or `$HOME/.local/bin`.
So if you want to use system commands like `pdflatex` or `x13as` you should either call them by their full path (e.g `/Library/TeX/texbin/pdflatex`)
or append the PATH in MATLAB by running `setenv('PATH', [getenv('PATH') ':/usr/local/bin:$HOME/.local/bin:/Library/TeX/texbin']);`.
Alternatively, you can create a `startup.m` file or change the system default PATH in the `/etc/paths` file.

### Optional: pass the full PATH to MATLAB to run system commands
If you start MATLAB from a terminal, you will get the PATH inherited from the shell.
However, when you click on the application icon in macOS, you are not running at the terminal level:
the program is run by launcher, which does not go through a shell login session.
In other words, you get the system default PATH which includes `/usr/bin:/bin:/usr/sbin:/sbin`, but not `/usr/local/bin` or `$HOME/.local/bin`.
So if you want to use system commands like `pdflatex` or `x13as` you should either call them by their full path (e.g `/Library/TeX/texbin/pdflatex`) or append the PATH in MATLAB by running `setenv('PATH', [getenv('PATH') ':/usr/local/bin:$HOME/.local/bin:/Library/TeX/texbin']);`.
Alternatively, you can create a `startup.m` file or change the system default PATH in the `/etc/paths` file.

Last tested on:
- macOS Ventura 13.3.1 (MacBook Air M1, MacBook Pro M2 MAX, M2 Virtual Machine using Parallels, Intel Virtual Machine using Quickemu)