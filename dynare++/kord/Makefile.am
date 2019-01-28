noinst_LIBRARIES = libkord.a

libkord_a_SOURCES = \
	approximation.cc \
	approximation.hh \
	decision_rule.cc \
	decision_rule.hh \
	dynamic_model.cc \
	dynamic_model.hh \
	faa_di_bruno.cc \
	faa_di_bruno.hh \
	first_order.cc \
	first_order.hh \
	global_check.cc \
	global_check.hh \
	kord_exception.hh \
	korder.cc \
	korder.hh \
	korder_stoch.cc \
	korder_stoch.hh \
	journal.cc \
	journal.hh \
	mersenne_twister.hh \
	normal_conjugate.cc \
	normal_conjugate.hh \
	random.cc \
	random.hh

libkord_a_CPPFLAGS = -I../sylv/cc -I../tl/cc -I../integ/cc -I$(top_srcdir)/mex/sources $(CPPFLAGS_MATIO)
libkord_a_CXXFLAGS = $(AM_CXXFLAGS) $(THREAD_CXXFLAGS)

check_PROGRAMS = tests

tests_SOURCES = tests.cc
tests_CPPFLAGS = -I../sylv/cc -I../tl/cc -I../integ/cc -I$(top_srcdir)/mex/sources
tests_CXXFLAGS = $(AM_CXXFLAGS) $(THREAD_CXXFLAGS)
tests_LDFLAGS = $(AM_LDFLAGS) $(LDFLAGS_MATIO)
tests_LDADD = libkord.a ../tl/cc/libtl.a ../sylv/cc/libsylv.a $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(LIBADD_MATIO)

check-local:
	./tests

CLEANFILES = out.txt
