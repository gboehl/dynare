RINTERNALS=/usr/share/R/include/

sylvcppsource := $(wildcard ../../sylv/cc/*.cc)
sylvhsource := $(wildcard ../../sylv/cc/*.hh)
sylvobjects := $(patsubst %.cc, %.o, $(sylvcppsource))

tlcppsource := $(wildcard ../../tl/cc/*.cc)
tlhsource := $(wildcard ../../tl/cc/*.hh)

kordcppsource := $(wildcard ../../kord/*.cc)
kordhsource := $(wildcard ../../kord/*.hh)

integcppsource := $(wildcard ../../integ/cc/*.cc)
integhsource := $(wildcard ../../integ/cc/*.hh)

parserhsource := $(wildcard ../../parser/cc/*.hh)
parsercppsource := $(wildcard ../parser/cc/*.cc)

utilshsource := $(wildcard ../../utils/cc/*.hh)
utilscppsource := $(wildcard ../utils/cc/*.cc)

srccpp := dynare3.cc dynare_model.cc planner_builder.cc dynare_atoms.cc dynare_params.cc  nlsolve.cc
objects := $(patsubst %.cc,../../src/%.o,$(srccpp)) \
$(patsubst %.yy,%_ll.o,$(wildcard ../../src/*.yy)) \
$(patsubst %.ll,%_tab.o,$(wildcard ../../src/*.ll))

PKG_CPPFLAGS= -I../../tl/cc -I../../sylv/cc -I../../kord -I../../src -I../.. -I$(RINTERNALS)
PKG_LIBS= ${LAPACK_LIBS} ${BLAS_LIBS} $(objects) $(kordobjects) $(integobjects) $(tlobjects) ../../parser/cc/parser.a ../../utils/cc/utils.a $(sylvobjects) -lpthread -llapack -lcblas -lf77blas -latlas -lg2c -lstdc++

ifneq ($(LD_LIBRARY_PATH),)	# use LD_LIBRARY_PATH from environment
	PKG_LIBS := -Wl,--library-path $(LD_LIBRARY_PATH) $(PKG_LIBS)
endif

dynareR.so: dynareR.o
	g++ -shared  -o dynareR.so dynareR.o -L/usr/lib/R/lib -lR $(PKG_LIBS)

dynareR.o: dynareR.cc
	g++ -I/usr/share/R/include -I/usr/share/R/include $(PKG_CPPFLAGS) \
	-fpic  -g -O2 -c dynareR.cc -o dynareR.o -DDEBUG

test: test.cc dynareR.cc
	g++ -O0 -g -o test test.cc -DDEBUG $(PKG_LIBS) $(PKG_CPPFLAGS)

test-debug:
	valgrind --leak-check=yes ./test
