check_PROGRAMS = tests

tests_SOURCES = tests.cc
tests_CPPFLAGS = -I../sylv/cc -I../tl/cc -I../integ/cc -I../utils/cc -I$(top_srcdir)/mex/sources
tests_LDFLAGS = $(AM_LDFLAGS) $(LDFLAGS_MATIO)
tests_LDADD = libkord.a ../tl/cc/libtl.a ../sylv/cc/libsylv.a ../utils/cc/libutils.a $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(LIBADD_MATIO)

check-local:
	./tests

CLEANFILES = out.txt
