# Configuration file for building Dynare Windows Installer
# Uses "NullSoft Scriptable Installer System", aka NSIS (see http://nsis.sourceforge.net)
# NSIS can be run from both Windows and Linux (see "nsis" package in Debian)

SetCompressor /SOLID lzma

Name "Dynare ${VERSION}"

OutFile "dynare-${VERSION}-win.exe"

InstallDir "c:\dynare\${VERSION}"

# Use the Modern User Interface (version 2)
!include MUI2.nsh

!define MUI_WELCOMEPAGE_TEXT "This wizard will guide you through the installation of Dynare ${VERSION}.$\n$\nDynare is distributed under the GNU General Public License (GPL) version 3.$\n$\nIf you accept the license, click Next button to continue the installation."
!insertmacro MUI_PAGE_WELCOME
!define MUI_COMPONENTSPAGE_NODESC
!insertmacro MUI_PAGE_COMPONENTS
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!define MUI_FINISHPAGE_LINK_LOCATION https://www.dynare.org
!define MUI_FINISHPAGE_LINK "Go to Dynare homepage"
!define MUI_FINISHPAGE_SHOWREADME $INSTDIR\README.txt
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_UNPAGE_WELCOME
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_UNPAGE_FINISH

!insertmacro MUI_LANGUAGE "English"

!define REGLOC "Software\Microsoft\Windows\CurrentVersion\Uninstall\Dynare ${VERSION}"
!define SMLOC "$SMPROGRAMS\Dynare ${VERSION}"

!macro DETERMINE_CONTEXT
 # Determine if we are admin or not
 # This will change the start menu directory and the registry root key (HKLM or HKLU)
 UserInfo::getAccountType
 Pop $0
 StrCmp $0 "Admin" +3
 SetShellVarContext current
 Goto +2
 SetShellVarContext all
!macroend


Section "Dynare core (preprocessor and M-files)"
 SectionIn RO
!insertmacro DETERMINE_CONTEXT
 SetOutPath $INSTDIR
 File README.txt ..\NEWS.md ..\license.txt

 SetOutPath $INSTDIR\matlab
 File /r ..\matlab\*.m ..\matlab\*.json ..\build-win-matlab\dynare_version.m

 SetOutPath $INSTDIR\preprocessor
 File ..\build-win-matlab\preprocessor\src\dynare-preprocessor.exe

 SetOutPath $INSTDIR\matlab\preprocessor64
 File ..\matlab\preprocessor64\dynare_m.exe

 SetOutPath $INSTDIR\matlab\dseries\externals\x13\windows\64
 File deps\lib64\x13as\x13as.exe

 SetOutPath $INSTDIR\contrib
 File /r ..\contrib\*.m

 SetOutPath $INSTDIR\scripts
 File ..\scripts\dynare.el

 WriteUninstaller $INSTDIR\uninstall.exe

 # Create entry in "Add/Remove programs"
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "DisplayName" "Dynare ${VERSION}"
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "DisplayVersion" "${VERSION}"
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "InstallLocation" $INSTDIR
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "UninstallString" "$INSTDIR\uninstall.exe"
 WriteRegDWORD SHELL_CONTEXT "${REGLOC}" "NoModify" 1
 WriteRegDWORD SHELL_CONTEXT "${REGLOC}" "NoRepair" 1

 # Remove "Modify" permission to "Authenticated Users".
 # This permission is inherited by default, since we're writing to "c:\dynare\";
 # it implies that anybody on the system can modify Dynare installation files,
 # which is a security issue.
 # Intigriti reference: DYNARE-7TULYVR8
 # On icacls, see e.g.: https://4sysops.com/archives/icacls-list-set-grant-remove-and-deny-permissions/#inheriting-permissions
 nsExec::Exec 'icacls "$INSTDIR" /inheritance:d'
 nsExec::Exec 'icacls "$INSTDIR" /remove:g "NT AUTHORITY\Authenticated Users"'
SectionEnd


Section "MEX files for MATLAB 64-bit, version 9.5 to 23.2 (R2018b to R2024a)"
 SetOutPath $INSTDIR\mex\matlab\win64-9.5-24.1
 File ..\build-win-matlab\*.mexw64
SectionEnd

Section "MinGW compiler for MATLAB 64-bit"
 SetOutPath $INSTDIR\mingw64
 File /r deps\mingw64\*
SectionEnd

Section "MEX files for Octave 9.1.0 (64-bit)"
 SetOutPath $INSTDIR\mex\octave\win64
 File ..\build-win-octave\*.mex
SectionEnd

Section "Documentation and examples"
 SetOutPath $INSTDIR\doc
 File ..\build-doc\dynare-manual.pdf ..\build-doc\guide.pdf ..\build-doc\bvar-a-la-sims.pdf ..\build-doc\dr.pdf ..\build-doc\sylvester.pdf ..\build-doc\tl.pdf ..\build-doc\preprocessor\doc\macroprocessor.pdf ..\build-doc\preprocessor\doc\preprocessor.pdf ..\build-doc\parallel.pdf ..\build-doc\gsa.pdf ..\build-doc\dseriesReporting.pdf

 SetOutPath $INSTDIR\doc\dynare-manual.html
 File /r ..\build-doc\dynare-manual.html\*

 SetOutPath $INSTDIR\examples
 File ..\examples\*.mod ..\examples\*.m

 # Create start menu entries
 CreateDirectory "${SMLOC}"
 CreateShortcut "${SMLOC}\Documentation.lnk" "$INSTDIR\doc"
 CreateShortcut "${SMLOC}\Examples.lnk" "$INSTDIR\examples"
SectionEnd


Section "Uninstall"
!insertmacro DETERMINE_CONTEXT

 # First delete the uninstaller
 Delete $INSTDIR\uninstall.exe
 Delete $INSTDIR\README.txt
 Delete $INSTDIR\NEWS.md
 Delete $INSTDIR\license.txt
 Delete $INSTDIR\VERSION
 Rmdir /r $INSTDIR\matlab
 Rmdir /r $INSTDIR\preprocessor
 Rmdir /r $INSTDIR\contrib
 Rmdir /r $INSTDIR\mex
 Rmdir /r $INSTDIR\doc
 Rmdir /r $INSTDIR\examples
 Rmdir /r $INSTDIR\scripts
 Rmdir /r $INSTDIR\mingw64
 # We don't force deletion of installation directory (with /r), to avoid deleting important files
 Rmdir $INSTDIR

 # Delete start menu entries
 Rmdir /r "${SMLOC}"

 # Delete entry in "Add/Remove programs"
 DeleteRegKey SHELL_CONTEXT "${REGLOC}"
SectionEnd
