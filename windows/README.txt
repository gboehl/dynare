Dynare
━━━━━━

For information about how to use Dynare, you should have a look at the
documentation located in the ‘doc’ subdirectory of your Dynare installation (you
should have a shortcut in your Start Menu to access it directly).

Beginners should start with the tutorial (under PDF format in ‘guide.pdf’).
There is also a complete reference manual documenting all Dynare functions
(under HTML format in ‘dynare-manual.html\index.html’, under PDF format in
‘dynare-manual.pdf’).

You can also get more information on the Dynare homepage:

  https://www.dynare.org

Note: Dynare comes with an automated uninstaller, which you can run from the
“Add or remove programs” menu of the Control Panel.


Using Dynare with MATLAB®
─────────────────────────

Dynare works on top of MATLAB®, any version ranging from 9.5 (R2018b) to 24.1
(R2024a). Only 64-bit versions are supported.

To use Dynare, you just have to add the ‘matlab’ subdirectory of your Dynare
installation to MATLAB® path. You have two options for doing that:

— Use the addpath command, by typing the following (assuming that you have
  installed Dynare at the standard location, and replacing ‘x.y’ by the correct
  version number):

    addpath c:\dynare\x.y\matlab

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

This version of Dynare is compiled for Octave 9.1.0 (MinGW, 64-bit),
and may not work with other versions of Octave. The recommended version of
Octave can be downloaded at:
  https://ftpmirror.gnu.org/gnu/octave/windows/octave-9.1.0-w64-installer.exe

Every time you run Octave, you should type the following command (assuming that
you have installed Dynare at the standard location, and replacing ‘x.y’ by
the correct version number):

  addpath c:\dynare\x.y\matlab

Note: if you don't want to type this command every time you run Octave, you can
put it in a file called ‘.octaverc’ in your home directory
(typically ‘c:\Users\USERNAME’). This file is run by Octave at every startup.

You can test your installation by typing ‘dynare’ at the Octave prompt. This
should give you an error message complaining that you did not specify a MOD
file.


Credits
───────

MATLAB® is a registered trademark of The Mathworks, Inc.
