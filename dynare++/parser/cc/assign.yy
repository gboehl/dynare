// -*- C++ -*-
/* Copyright © 2006-2011, Ondra Kamenik */

%code requires
{
#include "location.hh"
#define ASGN_LTYPE ogp::location_type
}

%code
{
#include "atom_assignings.hh"
#include <string>

void asgn_error(std::string);
int asgn_lex();
extern ogp::AtomAssignings* aparser;
}

%union
{
  int integer;
  char *string;
  char character;
}

%token EQUAL_SIGN SEMICOLON CHARACTER BLANK
%token <string> NAME;

%define api.prefix {asgn_}

%locations
%defines
%define parse.error verbose

%%

root : assignments | %empty;

assignments : assignments BLANK | assignments assignment | assignment | BLANK;

assignment : NAME EQUAL_SIGN material SEMICOLON {
	aparser->add_assignment(@1.off, $1, @1.ll, @3.off-@1.off, @3.ll + @4.ll);}
  | NAME space EQUAL_SIGN material SEMICOLON {
	aparser->add_assignment(@1.off, $1, @1.ll, @4.off-@1.off, @4.ll + @5.ll);}
  ;

material : material CHARACTER | material NAME | material BLANK | NAME | CHARACTER | BLANK;

space : space BLANK | BLANK;

%%

void
asgn_error(std::string mes)
{
  aparser->error(std::move(mes));
}
