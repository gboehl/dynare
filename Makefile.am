SUBDIRS = preprocessor doc tests mex/sources

if ENABLE_DYNAREPLUSPLUS
SUBDIRS += dynare++
endif

if ENABLE_MATLAB
SUBDIRS += mex/build/matlab
endif
if ENABLE_OCTAVE
SUBDIRS += mex/build/octave
endif

ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = \
	matlab \
	contrib \
	NEWS.md \
	license.txt \
	README.md \
	COPYING \
	CONTRIBUTING.md \
	windows \
	macOS \
	examples \
	scripts \
	.dir-locals.el

all-local:
# Create top-level preprocessor symlink needed by matlab/dynare.m (when Dynare is run
# from the locally-built copy)
	$(LN_S) -f $(abs_builddir)/preprocessor/src/dynare-preprocessor$(EXEEXT) $(abs_builddir)/preprocessor/dynare-preprocessor$(EXEEXT)
# Create backward-compatibility symlink for old location of preprocessor
	$(MKDIR_P) matlab/preprocessor64/
	$(LN_S) -f $(abs_builddir)/preprocessor/src/dynare-preprocessor$(EXEEXT) $(abs_srcdir)/matlab/preprocessor64/dynare_m$(EXEEXT)

clean-local:
	rm -f preprocessor/dynare-preprocessor$(EXEEXT)
	rm -rf matlab/preprocessor64/

dist-hook:
	rm -rf `find $(distdir)/matlab $(distdir)/examples -name *~`
	rm -f $(distdir)/matlab/dynare_version.m
	rm -rf $(distdir)/matlab/preprocessor64/
	$(MKDIR_P) $(distdir)/mex/matlab $(distdir)/mex/octave
	rm -rf `find $(distdir) -name '.git*'`

install-exec-local:
	$(MKDIR_P) $(DESTDIR)$(pkglibdir)/preprocessor
# The following will break under Windows, but we don’t use the install rule there
	$(LN_S) -f $(bindir)/dynare-preprocessor$(EXEEXT) $(DESTDIR)$(pkglibdir)/preprocessor/dynare-preprocessor$(EXEEXT)
	$(MKDIR_P) $(DESTDIR)$(pkglibdir)/contrib/ms-sbvar/TZcode
	cp -r contrib/ms-sbvar/TZcode/MatlabFiles $(DESTDIR)$(pkglibdir)/contrib/ms-sbvar/TZcode
	$(MKDIR_P) $(DESTDIR)$(pkglibdir)/contrib/jsonlab
	cp -r contrib/jsonlab $(DESTDIR)$(pkglibdir)/contrib/jsonlab
	cp -r examples $(DESTDIR)$(pkglibdir)
	cp -r matlab $(DESTDIR)$(pkglibdir)
	find $(DESTDIR)$(pkglibdir) -name LICENSE.md -delete
# Recreate backward-compatibility symlink
	rm -f $(DESTDIR)$(pkglibdir)/matlab/preprocessor64/dynare_m$(EXEEXT)
	$(LN_S) -f $(bindir)/dynare-preprocessor$(EXEEXT) $(DESTDIR)$(pkglibdir)/matlab/preprocessor64/dynare_m$(EXEEXT)
# We don’t fail over doc install rules, since the user may deliberately not have compiled them
# (e.g. with the “nodoc” option under Debian)
	$(MKDIR_P) $(DESTDIR)$(docdir)
	-cp doc/*.pdf doc/gsa/gsa.pdf doc/parallel/parallel.pdf doc/dseries-and-reporting/dseriesReporting.pdf preprocessor/doc/preprocessor/preprocessor.pdf preprocessor/doc/macroprocessor/macroprocessor.pdf doc/manual/build/latex/dynare-manual.pdf $(DESTDIR)$(docdir)
	-cp -r doc/manual/build/html $(DESTDIR)$(docdir)/dynare-manual.html
	$(MKDIR_P) $(DESTDIR)$(docdir)/dynare++
	-cp dynare++/doc/*.pdf $(DESTDIR)$(docdir)/dynare++

uninstall-local:
	rm -rf $(DESTDIR)$(pkglibdir)
	rm -rf $(DESTDIR)$(docdir)
