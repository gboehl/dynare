noinst_LIBRARIES = libinteg.a

libinteg_a_SOURCES = \
	quadrature.cc \
	quadrature.hh \
	quasi_mcarlo.cc \
	quasi_mcarlo.hh \
	product.cc \
	product.hh \
	smolyak.cc \
	smolyak.hh \
	vector_function.cc \
	vector_function.hh \
	precalc_quadrature.dat
libinteg_a_CPPFLAGS = -I../../sylv/cc -I../../tl/cc -I$(top_srcdir)/mex/sources
libinteg_a_CXXFLAGS = $(AM_CXXFLAGS) $(PTHREAD_CFLAGS)
