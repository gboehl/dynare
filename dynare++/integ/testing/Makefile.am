check_PROGRAMS = tests

tests_SOURCES = tests.cc
tests_CPPFLAGS = -I../cc -I../../tl/cc -I../../sylv/cc -I../../utils/cc -I$(top_srcdir)/mex/sources
tests_CXXFLAGS = $(AM_CXXFLAGS) $(THREAD_CXXFLAGS)
tests_LDFLAGS = $(AM_LDFLAGS) $(LDFLAGS_MATIO)
tests_LDADD = ../../tl/cc/libtl.a ../../sylv/cc/libsylv.a ../cc/libinteg.a ../../utils/cc/libutils.a $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(LIBADD_MATIO)

check-local:
	./tests
