noinst_LIBRARIES = libparser.a

GENERATED_FILES = assign_tab.cc formula_tab.cc matrix_tab.cc assign_tab.hh formula_tab.hh matrix_tab.hh assign_ll.cc formula_ll.cc matrix_ll.cc

libparser_a_SOURCES = \
	location.hh \
	atom_assignings.cc \
	atom_assignings.hh \
	atom_substitutions.cc \
	atom_substitutions.hh \
	dynamic_atoms.cc \
	dynamic_atoms.hh \
	fine_atoms.cc \
	fine_atoms.hh \
	formula_parser.cc \
	formula_parser.hh \
	matrix_parser.cc \
	matrix_parser.hh \
	parser_exception.cc \
	parser_exception.hh \
	static_atoms.cc \
	static_atoms.hh \
	static_fine_atoms.cc \
	static_fine_atoms.hh \
	tree.cc \
	tree.hh \
	$(GENERATED_FILES)

libparser_a_CPPFLAGS = -I../.. $(BOOST_CPPFLAGS)

BUILT_SOURCES = $(GENERATED_FILES)

EXTRA_DIST = assign.yy formula.yy matrix.yy assign.ll formula.ll matrix.ll

%_tab.cc %_tab.hh: %.yy
	$(YACC) -W -o$*_tab.cc $<

%_tab.$(OBJEXT): CXXFLAGS += -Wno-old-style-cast

%_ll.cc: %.ll
	$(LEX) -i -o$@ $<

%_ll.$(OBJEXT): CXXFLAGS += -Wno-old-style-cast
