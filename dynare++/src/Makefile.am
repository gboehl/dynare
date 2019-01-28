bin_PROGRAMS = dynare++

GENERATED_FILES = dynglob_ll.cc dynglob_tab.cc dynglob_tab.hh

dynare___SOURCES = \
	main.cc \
	dynare3.cc \
	dynare_atoms.hh \
	dynare_model.hh \
	forw_subst_builder.hh \
	planner_builder.cc \
	dynare3.hh \
	dynare_exception.hh \
	dynare_params.cc \
	planner_builder.hh \
	dynare_atoms.cc \
	dynare_model.cc \
	dynare_params.hh \
	forw_subst_builder.cc \
	nlsolve.cc \
	nlsolve.hh \
	$(GENERATED_FILES)

dynare___CPPFLAGS = -I../sylv/cc -I../tl/cc -I../kord -I../integ/cc -I.. -I$(top_srcdir)/mex/sources -DDYNVERSION=\"$(PACKAGE_VERSION)\" $(BOOST_CPPFLAGS) $(CPPFLAGS_MATIO)
dynare___LDFLAGS = $(AM_LDFLAGS) $(LDFLAGS_MATIO) $(BOOST_LDFLAGS)
dynare___LDADD = ../kord/libkord.a ../integ/cc/libinteg.a ../tl/cc/libtl.a ../parser/cc/libparser.a ../utils/cc/libutils.a ../sylv/cc/libsylv.a $(LIBADD_MATIO) $(noinst_LIBRARIES) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
dynare___CXXFLAGS = $(AM_CXXFLAGS) $(THREAD_CXXFLAGS)

BUILT_SOURCES = $(GENERATED_FILES)
EXTRA_DIST = dynglob.ll dynglob.yy

dynglob_tab.cc dynglob_tab.hh: dynglob.yy
	$(YACC) -d -odynglob_tab.cc dynglob.yy

dynglob_ll.cc: dynglob.ll
	$(LEX) -i -odynglob_ll.cc dynglob.ll
