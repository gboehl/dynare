check_PROGRAMS = tests

tests_SOURCES = factory.cc factory.hh monoms.cc monoms.hh tests.cc
tests_CPPFLAGS = -I../cc -I../../sylv/cc
tests_CXXFLAGS = $(AM_CXXFLAGS) $(PTHREAD_CFLAGS)
tests_LDFLAGS = $(AM_LDFLAGS) $(LDFLAGS_MATIO)
tests_LDADD = ../cc/libtl.a ../../sylv/cc/libsylv.a $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(PTHREAD_LIBS) $(LIBADD_MATIO)

check-local:
	./tests
