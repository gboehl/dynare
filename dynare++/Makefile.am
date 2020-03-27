SUBDIRS = utils/cc sylv parser/cc tl doc integ kord src tests

EXTRA_DIST = dynare_simul

install-exec-local:
	$(MKDIR_P) $(DESTDIR)$(pkglibdir)/dynare++
	cp -r dynare_simul/* $(DESTDIR)$(pkglibdir)/dynare++

uninstall-local:
	rm -rf $(DESTDIR)$(pkglibdir)/dynare++

