check_PROGRAMS = tests

tests_SOURCES = MMMatrix.cc MMMatrix.hh tests.cc
tests_LDADD = ../cc/libsylv.a ../../utils/cc/libutils.a $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
tests_CPPFLAGS = -I../cc -I../../utils/cc

EXTRA_DIST = *.mm

check-local:
	./tests
