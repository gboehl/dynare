// -*- C++ -*-
%code requires
{
#include "location.hh"
#define CSV_LTYPE ogp::location_type
}

%code
{
#include "csv_parser.hh"

void csv_error(const char*);
int csv_lex();
extern ogp::CSVParser* csv_parser;
}

%union
{
  char* string;
  int integer;
}

%token COMMA NEWLINE BOGUS
%token <string> ITEM

%define api.prefix {csv_};

%locations
%defines
%define parse.error verbose

%%

csv_file : line_list | line_list line;

line_list : line_list line newline | line newline | line_list newline | newline;

line : line comma | line item | item | comma;

comma : COMMA {csv_parser->nextcol();};

newline : NEWLINE {csv_parser->nextrow();};

item : ITEM {csv_parser->item(@1.off, @1.ll);};


%%

void
csv_error(const char* mes)
{
  csv_parser->csv_error(mes);
}
