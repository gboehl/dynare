.. default-domain:: dynare

##############################
Installation and configuration
##############################

Software requirements
=====================

Packaged versions of Dynare are available for Windows XP/Vista/7/8/10,
`Debian GNU/Linux <http://www.debian.org/>`_, `Ubuntu`_ and macOS 10.8
or later. Dynare should work on other systems, but some compilation
steps are necessary in that case.

In order to run Dynare, you need one of the following:

* MATLAB version 7.5 (R2007b) or above (MATLAB R2009b 64-bit for macOS);
* GNU Octave version 4.2.1 or above.

Packages of GNU Octave can be downloaded on the `Dynare website`_.

The following optional extensions are also useful to benefit from
extra features, but are in no way required:

* If under MATLAB: the Optimization Toolbox, the Statistics Toolbox,
  the Control System Toolbox;

* If under GNU Octave, the following `Octave-Forge
  <http://octave.sourceforge.net/>`_ packages: ``optim, io,
  statistics, control``.


Installation of Dynare
======================

After installation, Dynare can be used in any directory on your
computer. It is best practice to keep your model files in directories
different from the one containing the Dynare toolbox. That way you can
upgrade Dynare and discard the previous version without having to
worry about your own files.


On Windows
----------

Execute the automated installer called ``dynare-4.x.y-win.exe`` (where
``4.x.y`` is the version number), and follow the instructions. The
default installation directory is ``c:\dynare\4.x.y``.

After installation, this directory will contain several
sub-directories, among which are ``matlab``, ``mex`` and ``doc``.

The installer will also add an entry in your Start Menu with a
shortcut to the documentation files and uninstaller.

Note that you can have several versions of Dynare coexisting (for
example in ``c:\dynare``), as long as you correctly adjust your path
settings (see see :ref:`words-warning`).


On Debian GNU/Linux and Ubuntu
------------------------------

Please refer to the `Dynare wiki`_ for detailed instructions.

Dynare will be installed under ``/usr/lib/dynare``. Documentation will
be under ``/usr/share/doc/dynare-doc``.


On macOS
--------

To install Dynare for use with Matlab, execute the automated installer
called ``dynare-4.x.y.pkg`` (where *4.x.y* is the version number), and
follow the instructions. The default installation directory is
``/Applications/Dynare/4.x.y`` (please refer to the `Dynare wiki`_ for
detailed instructions).

After installation, this directory will contain several
sub-directories, among which are ``matlab``, ``mex`` and ``doc``.

Note that several versions of Dynare can coexist (by default in
``/Applications/Dynare``), as long as you correctly adjust your path
settings (see :ref:`words-warning`).

To install Dynare for Octave, first install Homebrew following the
instructions on their site: `https://brew.sh/
<https://brew.sh/>`_. Then install Octave, issuing the command ``brew
install octave`` at the Terminal prompt. You can then install the
latest stable version of Dynare by typing ``brew install dynare`` at
the Terminal prompt. You can also pass options to the installation
command. These options can be viewed by typing ``brew info dynare`` at
the Terminal prompt.


For other systems
-----------------

You need to download Dynare source code from the `Dynare website`_ and
unpack it somewhere.

Then you will need to recompile the pre-processor and the dynamic
loadable libraries. Please refer to `README.md
<https://git.dynare.org/Dynare/dynare/blob/master/README.md>`_.

.. _compil-install:

Compiler installation
=====================

Prerequisites on Windows
------------------------

If you are using MATLAB under Windows, install a C++ compiler on your
machine and configure it with MATLAB. There are at least two free
compilers you can use. First, there is Microsoft’s Visual Studio
Community (`https://www.visualstudio.com/`), which has the largest
history of MATLAB support, but requires much space on the
hard-disk. Second, since MATLAB R2015b, MATLAB supports the MinGW-w64
C/C++ Compiler from TDM-GCC. To install this compiler, use the Add-Ons
menu of MATLAB. Search for MinGW or select it from Features.

For older version of MATLAB, in particular before R2014a, it may
sometimes make sense to use the ``gcc`` compiler provided by
Cygwin. However, integrating it in MATLAB can be quite cumbersome and
should be considered as a legacy option. For details, see
`instructions on the Dynare wiki`_.

Prerequisites on Debian GNU/Linux and Ubuntu
--------------------------------------------

Users of MATLAB under Linux need to have a working compilation
environment installed. If not already present, it can be installed via
``apt-get install build-essential``.

Users of Octave under Linux should install the package for MEX file
compilation (under Debian or Ubuntu, it is called ``liboctave-dev``).

Prerequisites on macOS
----------------------

If you are using MATLAB under macOS, you should install the latest
version of XCode: `see instructions on the Dynare wiki
<https://git.dynare.org/Dynare/dynare/wikis/Install-on-MacOS>`_.


Configuration
=============

For MATLAB
----------

.. highlight:: matlab

You need to add the ``matlab`` subdirectory of your Dynare
installation to MATLAB path. You have two options for doing that:


* Using the ``addpath`` command in the MATLAB command window:

  Under Windows, assuming that you have installed Dynare in the
  standard location, and replacing ``4.x.y`` with the correct version
  number, type::

    >> addpath c:/dynare/4.x.y/matlab

  Under Debian GNU/Linux or Ubuntu, type::

    >> addpath /usr/lib/dynare/matlab

  Under macOS, assuming that you have installed Dynare in the standard
  location, and replacing ``4.x.y`` with the correct version number,
  type::

    >> addpath /Applications/Dynare/4.x.y/matlab

  MATLAB will not remember this setting next time you run it, and you
  will have to do it again.

* Via the menu entries:

  Select the “Set Path” entry in the “File” menu, then click on “Add
  Folder…”, and select the ``matlab`` subdirectory of ‘your Dynare
  installation. Note that you *should not* use “Add with
  Subfolders…”. Apply the settings by clicking on “Save”. Note that
  MATLAB will remember this setting next time you run it.


For GNU Octave
--------------

You need to add the ``matlab`` subdirectory of your Dynare
installation to Octave path, using the ``addpath`` at the Octave
command prompt.

Under Windows, assuming that you have installed Dynare in the standard
location, and replacing “*4.x.y*” with the correct version number,
type::

  octave:1> addpath c:/dynare/4.x.y/matlab

Under Debian GNU/Linux or Ubuntu, there is no need to use the
``addpath`` command; the packaging does it for you.

Under macOS, assuming that you have installed Dynare and Octave via
Homebrew, type::

  octave:1> addpath /usr/local/opt/dynare/lib/dynare/matlab

If you don’t want to type this command every time you run Octave, you
can put it in a file called ``.octaverc`` in your home directory
(under Windows this will generally be ``c:\Documents and
Settings\USERNAME\`` while under macOS it is
``/Users/USERNAME/``). This file is run by Octave at every startup.


.. _words-warning:

Some words of warning
---------------------

You should be very careful about the content of your MATLAB or Octave
path. You can display its content by simply typing ``path`` in the
command window.

The path should normally contain system directories of MATLAB or
Octave, and some subdirectories of your Dynare installation. You have
to manually add the ``matlab`` subdirectory, and Dynare will
automatically add a few other subdirectories at runtime (depending on
your configuration). You must verify that there is no directory coming
from another version of Dynare than the one you are planning to use.

You have to be aware that adding other directories (on top of the
dynare folders) to your MATLAB or Octave path can potentially create
problems if any of your M-files have the same name as a Dynare
file. Your routine would then override the Dynare routine, making
Dynare unusable.


.. warning::

   Never add all the subdirectories of the ``matlab`` folder to the
   MATLAB or Octave path. You must let Dynare decide which subdirectories
   have to be added to the MATLAB or Octave path. Otherwise, you may
   end up with a non optimal or un-usable installation of Dynare.


.. _Ubuntu: http://www.ubuntu.com/
.. _Dynare website: https://www.dynare.org/
.. _Dynare wiki: https://git.dynare.org/Dynare/dynare/wikis
.. _instructions on the Dynare wiki : http://www.dynare.org/DynareWiki/ConfigureMatlabWindowsForMexCompilation
