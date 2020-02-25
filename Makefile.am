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

all-local: preprocessor/src/dynare_m$(EXEEXT)
	if file preprocessor/src/dynare_m$(EXEEXT) | grep -q x86.64; then \
	  ARCH="64"; \
	else \
	  ARCH="32"; \
	fi; \
	mkdir -p $(abs_srcdir)/matlab/preprocessor$$ARCH && \
	$(LN_S) -f $(abs_builddir)/preprocessor/src/dynare_m$(EXEEXT) $(abs_srcdir)/matlab/preprocessor$$ARCH

clean-local:
	rm -rf $(abs_srcdir)/matlab/preprocessor32 $(abs_srcdir)/matlab/preprocessor64

dist-hook:
	rm -rf `find $(distdir)/matlab $(distdir)/examples -name *~`
	rm -rf $(distdir)/matlab/preprocessor* $(distdir)/matlab/dynare_version.m
	$(MKDIR_P) $(distdir)/mex/matlab $(distdir)/mex/octave
	rm -rf `find $(distdir) -name '.git*'`

install-exec-local:
	$(MKDIR_P) $(DESTDIR)$(pkglibdir)/contrib/ms-sbvar/TZcode
	cp -r contrib/ms-sbvar/TZcode/MatlabFiles $(DESTDIR)$(pkglibdir)/contrib/ms-sbvar/TZcode
	$(MKDIR_P) $(DESTDIR)$(pkglibdir)/contrib/jsonlab
	cp -r contrib/jsonlab $(DESTDIR)$(pkglibdir)/contrib/jsonlab
	cp -r examples $(DESTDIR)$(pkglibdir)
	cp -r matlab $(DESTDIR)$(pkglibdir)
	find $(DESTDIR)$(pkglibdir) -name LICENSE.md -delete
	rm -rf $(DESTDIR)$(pkglibdir)/matlab/preprocessor*
	if file preprocessor/src/dynare_m | grep -q x86.64; then \
	  ARCH="64"; \
	else \
	  ARCH="32"; \
	fi; \
	mkdir -p $(DESTDIR)$(pkglibdir)/matlab/preprocessor$$ARCH && \
	cp preprocessor/src/dynare_m $(DESTDIR)$(pkglibdir)/matlab/preprocessor$$ARCH

uninstall-local:
	rm -rf $(DESTDIR)$(pkglibdir)
