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
1. [**Docker**](#docker)

## General Instructions

### Prerequisites

A number of tools and libraries are needed in order to recompile everything. You don't necessarily need to install everything, depending on what you want to compile.

- The [GNU Compiler Collection](https://gcc.gnu.org/), version 10 or later, with
  gcc, g++ and gfortran
- [MATLAB](https://mathworks.com) (if you want to compile the MEX for MATLAB)
- [GNU Octave](https://www.octave.org) with
  - the development headers (if you want to compile the MEX for Octave)
  - the development libraries corresponding to the [UMFPACK](https://people.engr.tamu.edu/davis/suitesparse.html) packaged with Octave (if you want to compile the MEX for Octave)
  - the [statistics](https://octave.sourceforge.io/statistics/) package and, optionally, the [control](https://octave.sourceforge.io/control/), [io](https://octave.sourceforge.io/io/) and [optimization](https://octave.sourceforge.io/optim/) packages, either installed via your package manager or through [Octave Forge](https://octave.sourceforge.io/)
- [Meson](https://mesonbuild.com), version 0.64.0 or later
- [Pkgconf](http://pkgconf.org/), or another pkg-config implementation
- [Bash](https://www.gnu.org/software/bash/)
- [Boost libraries](https://www.boost.org), version 1.36 or later
- [Bison](https://www.gnu.org/software/bison/), version 3.2 or later
- [Flex](https://github.com/westes/flex), version 2.5.4 or later
- [MAT File I/O library](https://sourceforge.net/projects/matio/), version 1.5 or later (only when compiling for Octave)
- [SLICOT](http://www.slicot.org)
- [GSL library](https://www.gnu.org/software/gsl/)
- A decent LaTeX distribution (if you want to compile PDF documentation),
  ideally with Beamer
- For building the reference manual:
  - [Sphinx](https://www.sphinx-doc.org/)
  - [MathJax](https://www.mathjax.org/)
- [X-13ARIMA-SEATS Seasonal Adjustment Program](https://www.census.gov/data/software/x13as.html)

### Preparing the sources

If you have downloaded the sources from an official source archive or the source snapshot, just unpack it.

If you want to use Git, do the following from a terminal (note that you must
have the [Git LFS](https://git-lfs.github.com/) extension installed):
```sh
git clone --recurse-submodules https://git.dynare.org/Dynare/dynare.git
cd dynare
```
If you want a certain version (e.g. 5.x) , then add `--single-branch --branch 5.x` to the git clone command.

### Configuring the build directory

If you want to compile for MATLAB, please run the following (after adapting the path to MATLAB):
```sh
meson setup -Dmatlab_path=/usr/local/MATLAB/R2023b -Dbuildtype=debugoptimized build-matlab
```
The build directory will thus be `build-matlab`.

Or for Octave:
```sh
meson setup -Dbuild_for=octave -Dbuildtype=debugoptimized build-octave
```
The build directory will thus be `build-octave`.

Note that if you do not chose `build-matlab` (under MATLAB) or `build-octave`
(under Octave) as the build directory, you will need to set the environment
variable `DYNARE_BUILD_DIR` to the full path of your build tree, before running
MATLAB or Octave, if you want Dynare to be able to find the preprocessor and
the MEX files.

It is possible to specify various Meson options, see the Meson documentation
for more details. Modifying options of an existing build directory can be
done using the `meson configure` command.

### Building

For compiling the preprocessor and the MEX files:
```sh
meson compile -C <builddir>
```
where `<builddir>` is the build directory, typically either `build-matlab` or `build-octave`.

PDF and HTML documentation can be built with:
```sh
meson compile -C <builddir> doc
```

### Check

Dynare comes with unit tests (in the MATLAB functions) and integration tests (under the `tests` subfolder). All the tests can be run with:
```sh
meson test -C <builddir>
```

Depending on the performance of your machine, this can take several hours.

Note that running the testsuite with Octave requires the additional packages `pstoedit`, `epstool`, `xfig`, and `gnuplot`. 

Often, it does not make sense to run the complete testsuite. For instance, if you modify codes only related to the perfect foresight model solver, you can decide to run only a subset of the integration tests, with:
```sh
meson test -C <builddir> --suite deterministic_simulations
```
This will run all the integration tests in `tests/deterministic_simulations`.
This syntax also works with a nested directory (e.g. `--suite deterministic_simulations/purely_forward`).

Finally if you want to run a single integration test, e.g. `deterministic_simulations/lbj/rbc.mod`:
```sh
meson test -C <builddir> deterministic_simulations/lbj/rbc.mod
```
NB: Some individual tests cannot be run using that syntax, if they are a dependency in a chain of tests (see the `mod_and_m_tests` variable `meson.build`); in that case, you should use the name of the last `.mod` file in the chain as the test name to be passed to `meson test`.

## Debian or Ubuntu

All the prerequisites are packaged:

- `gcc`
- `g++`
- `gfortran`
- `octave-dev` (or `liboctave-dev` on older Debian/Ubuntu releases)
- `libboost-graph-dev`
- `libgsl-dev`
- `libmatio-dev`
- `libslicot-dev` and `libslicot-pic`
- `libsuitesparse-dev`
- `flex` and `libfl-dev`
- `bison`
- `meson`
- `pkgconf`
- `texlive`
- `texlive-publishers` (for Econometrica bibliographic style)
- `texlive-latex-extra` (for fullpage.sty)
- `texlive-fonts-extra` (for ccicons)
- `texlive-science` (for amstex)
- `lmodern` (for macroprocessor PDF)
- `python3-sphinx`
- `tex-gyre`
- `latexmk`
- `libjs-mathjax`
- `x13as`

You can install them all at once with:
```sh
apt install gcc g++ gfortran octave-dev libboost-graph-dev libgsl-dev libmatio-dev libslicot-dev libslicot-pic libsuitesparse-dev flex libfl-dev bison meson pkgconf texlive texlive-publishers texlive-latex-extra texlive-fonts-extra texlive-science lmodern python3-sphinx make tex-gyre latexmk libjs-mathjax x13as
```
If you use MATLAB, we strongly advise to also `apt install matlab-support` and confirm to rename the GCC libraries shipped with MATLAB to avoid possible conflicts with GCC libraries shipped by your distribution.

## Fedora, CentOS or RHEL

Almost all prerequisites are packaged:

- `gcc`, `gcc-c++`
- `gcc-gfortran`
- `boost-devel`
- `gsl-devel`
- `matio-devel`
- `suitesparse-devel`
- `flex`
- `bison`
- `meson`
- `redhat-rpm-config`
-  `octave`, `octave-devel`, `octave-statistics`, `octave-io`, `octave-optim`, `octave-control`
- `texlive-scheme-minimal`, `texlive-collection-publishers`, `texlive-collection-latexextra`, `texlive-collection-fontsextra`, `texlive-collection-latexrecommended`, `texlive-collection-science`, `texlive-collection-plaingeneric`, `texlive-lm`
- `python3-sphinx`
- `latexmk`
- `mathjax`
- `make` (for building Slicot)

You can install them all at once with:
```sh
# Minimal packages
dnf install -y gcc gcc-c++ make gcc-gfortran boost-devel gsl-devel matio-devel suitesparse-devel flex bison meson redhat-rpm-config
# Octave packages
dnf install octave octave-devel octave-statistics octave-io octave-optim octave-control
# Documentation packages (only needed if you build documentation)
dnf install texlive-scheme-minimal texlive-collection-publishers texlive-collection-latexextra texlive-collection-fontsextra texlive-collection-latexrecommended texlive-collection-science texlive-collection-plaingeneric texlive-lm python3-sphinx latexmk mathjax
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
The documentation packages have slightly different names in CentOS and RHEL, but this should only impact you if you build the documentation.

`Slicot` and `x13as` need to be compiled from source:

```sh
# compile slicot from source and put it into /home/$USER/dynare/slicot/lib/
mkdir -p /home/$USER/dynare/slicot
cd /home/$USER/dynare/slicot
wget https://deb.debian.org/debian/pool/main/s/slicot/slicot_5.0+20101122.orig.tar.gz
tar xf slicot_5.0+20101122.orig.tar.gz
cd slicot-5.0+20101122
mkdir -p /home/$USER/dynare/slicot/lib
# The following two lines are only for MATLAB
make FORTRAN=gfortran OPTS="-O2 -fPIC -fdefault-integer-8" LOADER=gfortran lib
cp slicot.a /home/$USER/dynare/slicot/lib/libslicot64_pic.a
# The following two lines are only for Octave
make FORTRAN=gfortran OPTS="-O2 -fPIC" LOADER=gfortran lib
cp slicot.a /home/$USER/dynare/slicot/lib/libslicot_pic.a

# compile x13as from source and put it into /usr/bin/
mkdir -p /home/$USER/dynare/x13as
cd /home/$USER/dynare/x13as
wget https://www2.census.gov/software/x-13arima-seats/x13as/unix-linux/program-archives/x13as_asciisrc-v1-1-b60.tar.gz
tar xf x13as_asciisrc-v1-1-b60.tar.gz
sed -i "s|-static| |" makefile.gf # this removes '-static' in the makefile.gf
make -f makefile.gf FFLAGS="-O2 -std=legacy" PROGRAM=x13as
sudo cp x13as /usr/bin/
```

If you use MATLAB, we strongly advise to also rename or exclude the GCC libraries shipped with MATLAB to avoid possible conflicts with GCC libraries shipped by Fedora, see e.g. [Matlab on Fedora 33](https://mutschler.eu/linux/install-guides/fedora-post-install/#matlab) or [MATLAB-ArchWiki](https://wiki.archlinux.org/index.php/MATLAB) for instructions.

Now use the following commands if using MATLAB (adapt them for Octave, see above):
```sh
cd /home/$USER/dynare
git clone --recurse-submodules https://git.dynare.org/dynare/dynare.git unstable
cd unstable
meson setup -Dmatlab_path=/usr/local/MATLAB/R2023b -Dfortran_args="[ '-B', '/home/$USER/dynare/slicot']" -Dbuildtype=debugoptimized build-matlab
meson compile -C build-matlab
```

If your distribution ships an older version of `bison`, compile it from source and append it *temporarily* to your path before running meson:
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
pacman -S git bison flex make tar mingw-w64-x86_64-meson mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-boost mingw-w64-x86_64-gsl mingw-w64-x86_64-matio mingw-w64-x86_64-pkgconf
```
- Compile and install SLICOT
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
```
- Configure Dynare from the source directory:
```sh
meson setup -Dmatlab_path=<…> -Dbuildtype=debugoptimized -Dprefer_static=true -Dfortran_args="['-B','/usr/local/lib']" build-matlab
```
where the path of MATLAB is specified. Note that you should use
the MSYS2 notation and not put spaces in the MATLAB path, so you probably want
to use something like `/c/Progra~1/MATLAB/…`. Alternatively, if your filesystem
does not have short filenames (8dot3), then you can run `mkdir -p
/usr/local/MATLAB && mount c:/Program\ Files/MATLAB /usr/local/MATLAB`, and
then pass `/usr/local/MATLAB/…` as MATLAB path to the configure script.
- Compile:
```sh
meson compile -C build-matlab
```
- Run the testsuite:
```sh
meson test -C build-matlab
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
create a temporary environment variable `BREWDIR` to ensure that the correct packages are used.

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
Create environment variables for which platform you want to compile for, i.e. either `arm64` or `x86_64`:

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

Install Homebrew using the environment variables:
```sh
arch -$ARCH /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
The prefix `arch -arm64` or `arch -x86_64` makes sure that you are installing the correct packages.
Don't forget to run the displayed commands (**Next steps**) in the terminal to add Homebrew to your PATH.

If you have both Homebrew installations installed, make sure that you are accessing the correct packages by temporarily (!) prepending it to the path:
```sh
export PATH="$BREWDIR/bin:$PATH"
```


- Install required Homebrew packages:
```sh
arch -$ARCH $BREWDIR/bin/brew install meson bison flex boost gcc gsl libmatio veclibfort octave sphinx-doc docutils wget pkg-config git-lfs
```
If you are installing `git-lfs` for the first time, you need to run `git lfs install` once after installing it.

- Link the sphinx-doc package to be able to compile the documentation:
```sh
arch -$ARCH $BREWDIR/bin/brew link --force sphinx-doc
```

- Install [MacTeX](http://www.tug.org/mactex/index.html) using the universal installer, if you want to build the documentation. MacTeX runs natively on both ARM and Intel machines. On Apple Silicon, it is advised to symlink `pdflatex`, `bibtex` and `latexmk` into `/usr/local/bin`:
```sh
sudo ln -s /Library/TeX/texbin/pdflatex /usr/local/bin/pdflatex
sudo ln -s /Library/TeX/texbin/bibtex /usr/local/bin/bibtex
sudo ln -s /Library/TeX/texbin/latexmk /usr/local/bin/latexmk
```
If you don't have admin privileges, then you can also symlink them into `$HOME/.local/bin` and add this folder to your PATH.

- Install MATLAB and additional toolboxes.
We recommend, but don't require, the following: Optimization, Global Optimization, Statistics and Machine Learning, Econometrics, and Control System.
For Apple Silicon: MATLAB offers a native Apple silicon version (arm64) as of version R2023b, see [the official instructions](https://de.mathworks.com/support/requirements/apple-silicon.html) how to install it.
You can also run the Intel version (x86_64) under Rosetta 2.
Don't forget to run MATLAB at least once to make sure you have a valid license.

- Create a folder for Dynare and its dependencies
```sh
export DYNAREDIR=$HOME/dynare
```

- Compile and install SLICOT
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
```

- Compile and install the X-13ARIMA-SEATS Seasonal Adjustment Program
```sh
mkdir -p $DYNAREDIR/x13as
cd $DYNAREDIR/x13as
curl -O https://www2.census.gov/software/x-13arima-seats/x13as/unix-linux/program-archives/x13as_asciisrc-v1-1-b60.tar.gz
tar xf x13as_asciisrc-v1-1-b60.tar.gz
sed -i '' 's/-static//g' makefile.gf
make -j$(sysctl -n hw.ncpu) -f makefile.gf FC=$BREWDIR/bin/gfortran LINKER=$BREWDIR/bin/gcc-13 FFLAGS="-O2 -std=legacy" LDFLAGS=-static-libgcc LIBS="$BREWDIR/lib/gcc/current/libgfortran.a /$BREWDIR/lib/gcc/current/libquadmath.a" PROGRAM=x13as
sudo cp $DYNAREDIR/x13as/x13as /usr/local/bin/x13as
cd $DYNAREDIR
x13as
```
Alternatively, if you don't have admin privileges you can install it into `$HOME/.local/bin` and add this folder to your PATH.

### Compile Dynare from source
The following commands will download the Dynare source code and compile
it. They should be entered at the command prompt in Terminal.app from the
folder where you want Dynare installed.

- Prepare the Dynare sources for the unstable version:
```sh
git clone --recurse-submodules https://git.dynare.org/Dynare/dynare.git $DYNAREDIR/unstable
cd $DYNAREDIR/unstable
```
If you want a certain version (e.g. 5.x) , then add `--single-branch --branch 5.x` to the git clone command.

- Configure Dynare from the source directory:
```sh
export BUILDDIR=build-matlab
export MATLABPATH=/Applications/MATLAB_R2023b.app
arch -$ARCH meson setup --native-file macOS/homebrew-native-$ARCH.ini -Dmatlab_path=$MATLABPATH -Dbuildtype=debugoptimized -Dfortran_args="['-B','$DYNAREDIR/slicot/lib']" $BUILDDIR
```
where you need to adapt the path to MATLAB.
Similarly, if you want to compile for Octave, replace the `-Dmatlab_path` option by `-Dbuild_for=octave`, and change the build directory to `build-octave`.

- Compile:
```sh
arch -$ARCH meson compile -C $BUILDDIR
```
If no errors occured, you are done. Dynare is now ready to use.

- If you additionally want to compile the documentation run:
```sh
arch -$ARCH meson compile -C $BUILDDIR doc
```

- Optionally, run the testsuite:
```sh
arch -$ARCH meson test -C $BUILDDIR --num-processes=$(sysctl -n hw.perflevel0.physicalcpu)
```
where `--num-processes` specifies the number of parallel processes to use for the testsuite (here set to the number of performance cores on your mac).

### Optional: pass the full PATH to MATLAB to run system commands
If you start MATLAB from a terminal, you will get the PATH inherited from the shell.
However, when you click on the application icon in macOS, you are not running at the terminal level:
the program is run by launcher, which does not go through a shell login session.
In other words, you get the system default PATH which includes `/usr/bin:/bin:/usr/sbin:/sbin`, but not `/usr/local/bin` or `$HOME/.local/bin`.
So if you want to use system commands like `pdflatex`, `latexmk` or `x13as` you should either call them by their full path (e.g `/Library/TeX/texbin/pdflatex`)
or append the PATH by running `setenv('PATH', [getenv('PATH') ':/usr/local/bin:$HOME/.local/bin:/Library/TeX/texbin']);` in your MATLAB command line once,
e.g. by adding this to your mod file. Alternatively, you can create a `startup.m` file or change the system default PATH in the `/etc/paths` file.


## Docker
We offer a variety of pre-configured Docker containers for Dynare, pre-configured with Octave and MATLAB including all recommended toolboxes.
These are readily available for your convenience on [Docker Hub](https://hub.docker.com/r/dynare/dynare).
The `scripts/docker` folder contains [information and instructions](scripts/docker/README.md) to interact, built and customize the containers.
