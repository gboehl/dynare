SUBDIRS = preprocessor doc tests mex/sources

if HAVE_BLAS
if HAVE_LAPACK
if HAVE_MATIO
SUBDIRS += dynare++
endif
endif
endif

# MEX must be built after dynare++ (because of kordepert)
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
	NEWS \
	license.txt \
	README.md \
	windows/dynare.nsi \
	windows/mexopts-win32.bat \
	windows/mexopts-win64.bat \
	windows/README.txt \
	osx \
	examples \
	scripts \
	.dir-locals.el

dist-hook:
	rm -rf `find $(distdir)/matlab $(distdir)/examples -name *~`
	rm -f $(distdir)/matlab/dynare_m$(EXEEXT) $(distdir)/matlab/dynare_version.m
	$(MKDIR_P) $(distdir)/mex/matlab $(distdir)/mex/octave
	rm -rf `find $(distdir)/contrib -name '.git*'`

install-exec-local:
	$(MKDIR_P) $(DESTDIR)$(pkglibdir)/contrib/ms-sbvar/TZcode
	cp -r examples $(DESTDIR)$(pkglibdir)
	cp -r matlab $(DESTDIR)$(pkglibdir)
	rm -f $(DESTDIR)$(pkglibdir)/matlab/dynare_m
	cp preprocessor/dynare_m $(DESTDIR)$(pkglibdir)/matlab
	cp -r contrib/ms-sbvar/TZcode/MatlabFiles $(DESTDIR)$(pkglibdir)/contrib/ms-sbvar/TZcode

uninstall-local:
	rm -f  $(DESTDIR)$(bindir)/dynare++
	rm -rf $(DESTDIR)$(pkglibdir)
