// -*- C++ -*-
// Copyright © 2007-2011, Ondra Kamenik

%code requires
{
#include "location.hh"
#define NAMELIST_LTYPE ogp::location_type
}

%code
{
#include "namelist.hh"

void namelist_error(const char*);
int namelist_lex();
extern ogp::NameListParser* name_list_parser;
}

%union
{
  int integer;
  char *string;
  char character;
}

%token COMMA CHARACTER
%token <string> NAME;

%define api.prefix {namelist_}

%locations
%defines
%define parse.error verbose

%%

namelist : namelist NAME       {name_list_parser->add_name($2);}
         | namelist COMMA NAME {name_list_parser->add_name($3);}
         | NAME                {name_list_parser->add_name($1);}
         ;

%%

void
namelist_error(const char* mes)
{
  name_list_parser->namelist_error(mes);
}
