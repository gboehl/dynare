Dynare
━━━━━━

For information about how to use Dynare, you should have a look at the
documentation located in the ‘doc’ subdirectory of your Dynare installation (you
should have a shortcut in your Start Menu to access it directly).

Beginners should start with the tutorial (under PDF format in ‘guide.pdf’).
There is also a complete reference manual documenting all Dynare functions
(under HTML format in ‘dynare-manual.html\index.html’, under PDF format in
‘dynare-manual.pdf’).

You can also get more information on the web, on Dynare homepage:

  https://www.dynare.org

Or on Dynare Wiki:

  https://www.dynare.org/DynareWiki

Note: Dynare comes with an automated uninstaller, which you can run from the
“Add/Remove Programs” menu of the Control Panel.


Using Dynare with MATLAB®
─────────────────────────

Dynare requires MATLAB® version 7.9 (R2009b) or above. With older versions of
MATLAB®, it may fail or give unexpected results.

To use Dynare, you just have to add the ‘matlab’ subdirectory of your Dynare
installation to MATLAB® path. You have two options for doing that:

— Use the addpath command, by typing the following (assuming that you have
  installed Dynare at the standard location, and replacing ‘4.x.y’ by the correct
  version number):

    addpath c:\dynare\4.x.y\matlab

  MATLAB® will not remember this setting next time you run it, and you will
  have to do it again.

— Select the “Set Path” entry in the “File” menu, then click on “Add Folder…”,
  and select the ‘matlab’ subdirectory of your Dynare installation. Note that
  you must not use “Add with Subfolders…”. Apply the settings by clicking on
  “Save”. Note that MATLAB® will remember this setting next time you run it.

You can test your installation by typing ‘dynare’ at the MATLAB® prompt. This
should give you an error message complaining that you did not specify a MOD
file.


Using Dynare with Octave
────────────────────────

Dynare also works on top of GNU Octave, a free clone of MATLAB® (see
<https://www.octave.org>).

This version of Dynare is compiled for Octave 4.4.1 (MinGW, 32-bit and 64-bit),
and may not work with other versions of Octave. The recommended version of
Octave can be downloaded at:

— For 64-bit systems:
  https://ftpmirror.gnu.org/gnu/octave/windows/octave-4.4.1-w64-installer.exe
— For 32-bit systems:
  https://ftpmirror.gnu.org/gnu/octave/windows/octave-4.4.1-w32-installer.exe

Every time you run Octave, you should type the two following commands (assuming
that you have installed Dynare at the standard location, and replacing ‘4.x.y’
by the correct version number):

  addpath c:\dynare\4.x.y\matlab

Note: if you don't want to type this command every time you run Octave, you can
put it in a file called ‘.octaverc’ in your home directory (‘c:\Documents and
Settings\USERNAME’ for Windows XP or ‘c:\Users\USERNAME’ for more recent
versions of Windows). This file is run by Octave at every startup.

You can test your installation by typing ‘dynare’ at the Octave prompt. This
should give you an error message complaining that you did not specify a MOD
file.


Dynamic Loadable Libraries
──────────────────────────

For better performance, some parts of Dynare are written in the C++ language,
which is faster than standard M-files. These parts are compiled and distributed
as dynamic loadable libraries (DLL), located in the 'mex' subdirectory of your
Dynare installation.

If the DLL are correctly detected by MATLAB® or Octave, the following should
be displayed when you launch Dynare:

  Configuring Dynare ...
  [mex] Generalized QZ.
  [mex] Sylvester equation solution.
  [mex] Kronecker products.
  [mex] Sparse kronecker products.
  [mex] Local state space iteration (second order).
  [mex] Bytecode evaluation.
  [mex] k-order perturbation solver.
  [mex] k-order solution simulation.
  [mex] Quasi Monte-Carlo sequence (Sobol).

On the contrary, if DLL are not detected, Dynare will fallback on
slower alternatives written in M-files (only for some of the DLLs),
and display the following:

  Configuring Dynare ...
  [m]   Generalized QZ.
  [m]   Sylvester equation solution.
  [m]   Kronecker products.
  [m]   Sparse kronecker products.
  [m]   Local state space iteration (second order).
  [no]  Bytecode evaluation.
  [no]  k-order perturbation solver.
  [no]  k-order solution simulation.
  [no]  Quasi Monte-Carlo sequence (Sobol).

In this last case, Dynare will run correctly on the basic features,
but with suboptimal speed, and some features will be missing. There
could be several reasons for MATLAB® or Octave failing to detect
the DLL:

— Your path settings may be wrong. Make sure that the ‘matlab’ subdirectory of
  your Dynare installation is the only Dynare directory present in the path
  variable.

— Your MATLAB® or Octave version may be incompatible with the provided
  binaries.

— You may have a custom M-file in your search path with the same name than a
  DLL, therefore overriding it.


Credits
───────

MATLAB® is a registered trademark of The Mathworks, Inc.
