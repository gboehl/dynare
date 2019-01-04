noinst_PROGRAMS = quadrature-points

quadrature_points_SOURCES = quadrature-points.cc
quadrature_points_CPPFLAGS = -I../.. -I../../sylv/cc -I../../integ/cc -I../../tl/cc
quadrature_points_CXXFLAGS = $(AM_CXXFLAGS) $(PTHREAD_CFLAGS)
quadrature_points_LDADD = ../cc/libinteg.a ../../tl/cc/libtl.a ../../parser/cc/libparser.a ../../sylv/cc/libsylv.a ../../utils/cc/libutils.a $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(PTHREAD_LIBS)
