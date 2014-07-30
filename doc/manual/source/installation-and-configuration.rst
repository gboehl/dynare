.. _installation-and-introduction:

********************************
 Installation and configuration
********************************

.. contents::
   :local:
   :depth: 2

.. installation_soft-requirements:

Software requirements
---------------------

Packaged versions of Dynare are available for Windows XP/Vista/7/8,
`Debian GNU/Linux`_, `Ubuntu`_ and Mac OS X Leopard/Snow Leopard. Dynare
should work on other systems, but some compilation steps are necessary
in that case.

In order to run Dynare, you need one of the following:

 * `MATLAB`_ version 7.5 (R2007b) or above;
 * `GNU Octave`_ version 3.6 or above.

Packages of GNU Octave can be downloaded on the `Dynare website`_.

The following optional extensions are also useful to benefit from
extra features, but are in no way required:

 * If under MATLAB: the optimization toolbox, the statistics toolbox,
   the control system toolbox;
 * If under GNU Octave, the following `Octave-Forge`_ packages: optim,
   io, java, statistics, control.

If you plan to use the ``use_dll`` option of the model command,
you will need to install the necessary requirements for compiling MEX
files on your machine. If you are using MATLAB under Windows, install
a C++ compiler on your machine and configure it with MATLAB: see
instructions on the `Dynare wiki`_. Users of Octave under Linux should
install the package for MEX file compilation (under Debian or Ubuntu,
it is called ‘liboctave-dev’). If you are using Octave or MATLAB under
Mac OS X, you should install the latest version of XCode: see
instructions on the Dynare wiki. Mac OS X Octave users will also need
to install gnuplot if they want graphing capabilities. Users of MATLAB
under Linux and Mac OS X, and users of Octave under Windows, normally
need to do nothing, since a working compilation environment is
available by default.


.. installation_install:

Installation of Dynare
----------------------

.. installation_install_windows:

On Windows
^^^^^^^^^^

Execute the automated installer called dynare-4.x.y-win.exe (where
4.x.y is the version number), and follow the instructions. The default
installation directory is c:\\dynare\\4.x.y. After installation, this
directory will contain several sub-directories, among which are
‘matlab’, ‘mex’ and ‘doc’. The installer will also add an entry in
your Start Menu with a shortcut to the documentation files and
uninstaller.

Note that you can have several versions of Dynare coexisting (for
example in ‘c:\\dynare’), as long as you correctly adjust your path
settings (see section Some words of warning).


.. installation_install_linux:

On Debian GNU/Linux and Ubuntu
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Please refer to the `Dynare wiki`_ for detailed instructions. Dynare
will be installed under ‘/usr/share/dynare’ and
‘/usr/lib/dynare’. Documentation will be under
‘/usr/share/doc/dynare’.


.. installation_install_macos:

On Mac OS X
^^^^^^^^^^^

Execute the automated installer called ‘dynare-4.x.y.pkg’ (where 4.x.y
is the version number), and follow the instructions. The default
installation directory is ‘/Applications/Dynare/4.x.y’.

Please refer to the `Dynare wiki`_ for detailed instructions. After
installation, this directory will contain several sub-directories,
among which are ‘matlab’, ‘mex’ and ‘doc’.

Note that you can have several versions of Dynare coexisting (for
example in ‘/Applications/Dynare’), as long as you correctly adjust
your path settings (see section Some words of warning).


.. installation_config:

Configuration
-------------


.. installation_config_matlab:

For MATLAB
^^^^^^^^^^

You need to add the ‘matlab’ subdirectory of your Dynare installation
to MATLAB path. You have two options for doing that:

* Using the addpath command in the MATLAB command window.

  * Under Windows, assuming that you have installed Dynare in the
    standard location, and replacing 4.x.y with the correct version
    number, type:::

    >> addpath c:\dynare\4.x.y\matlab

  * Under Debian GNU/Linux or Ubuntu, type:::

    >> addpath /usr/share/dynare/matlab

  * Under Mac OS X, assuming that you have installed Dynare in the
    standard location, and replacing 4.x.y with the correct version
    number, type:::1

    >> addpath /Applications/Dynare/4.x.y/matlab

 MATLAB will not remember this setting next time you run it, and you
 will have to do it again. The advantage of this option is that it
 ensures that Dynare is on top of the MATLAB path.

* Via the menu entries: Select the “Set Path” entry in the “File”
  menu, then click on “Add Folder…”, and select the ‘matlab’
  subdirectory of your Dynare installation. Note that you should not
  use “Add with Subfolders…”. Apply the settings by clicking on
  “Save”. Note that MATLAB will remember this setting next time you
  run it.


.. installation_config_octave:

For GNU Octave
^^^^^^^^^^^^^^

You need to add the ‘matlab’ subdirectory of your Dynare installation
to Octave path, using the addpath at the Octave command prompt.

* Under Windows, assuming that you have installed Dynare in the standard
  location, and replacing “4.x.y” with the correct version number, type:::

    >> addpath c:\dynare\4.x.y\matlab

* Under Debian GNU/Linux or Ubuntu, there is no need to use the
  addpath command; the packaging does it for you.

* Under Mac OS X, assuming that you have installed Dynare in the
  standard location, and replacing “4.x.y” with the correct version
  number, type:::

    >> addpath /Applications/Dynare/4.x.y/matlab

If you do not want to type this command every time you run Octave, you can put it in a file called ‘.octaverc’ in your home directory (under Windows this will generally be ‘c:\\Documents and Settings\\USERNAME\\’ while under Mac OS X it is ‘/Users/USERNAME/’). This file is run by Octave at every startup.


.. installation_config_warnings:

Some words of warning
^^^^^^^^^^^^^^^^^^^^^

You should be very careful about the content of your MATLAB or Octave
path. You can display its content by simply typing ``path`` in the command
window.

The path should normally contain system directories of MATLAB or
Octave, and some subdirectories of your Dynare installation. You have
to manually add the ‘matlab’ subdirectory, and Dynare will
automatically add a few other subdirectories at runtime (depending on
your configuration). You must verify that there is no directory coming
from another version of Dynare than the one you are planning to use.

You have to be aware that adding other directories to your path can
potentially create problems if any of your M-files have the same name
as a Dynare file. Your file would then override the Dynare file,
making Dynare unusable.

.. _Debian GNU/Linux: http://www.debian.org/
.. _Ubuntu: http://www.ubuntu.com/
.. _Dynare website: http://www.dynare.org/download/octave
.. _Octave-Forge: http://octave.sourceforge.net/

.. _GNU Octave: http://www.octave.org/
.. _MATLAB: http://www.mathworks.com/products/matlab/
.. _Dynare wiki: http://www.dynare.org/DynareWiki
