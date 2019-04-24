// -*- C++ -*-
// Copyright © 2006-2011, Ondra Kamenik

%code requires
{
#include "location.hh"
#define MATRIX_LTYPE ogp::location_type
}

%code
{
#include "matrix_parser.hh"

void matrix_error(std::string);
int matrix_lex();
extern ogp::MatrixParser* mparser;
}

%union
{
  double val;
  int integer;
}

%token NEW_ROW
%token <val> DNUMBER

%define api.prefix {matrix_};

%locations
%defines
%define parse.error verbose

%%

matrix : first_row other_rows
    | first_row other_rows empty_rows
    | first_row empty_rows other_rows empty_rows
    | first_row empty_rows other_rows
    | empty_rows first_row other_rows
    | empty_rows first_row other_rows empty_rows
    | empty_rows first_row empty_rows other_rows empty_rows
    | empty_rows first_row empty_rows
    | first_row empty_rows
    | empty_rows first_row
    | first_row
    | empty_rows
    ;

empty_rows : empty_rows NEW_ROW | NEW_ROW;

lod : DNUMBER {mparser->add_item($1);}
    | lod DNUMBER {mparser->add_item($2);}
    ;

first_row : lod;

other_rows : other_rows one_row | other_rows empty_rows one_row |one_row ;

one_row : NEW_ROW {mparser->start_row();} lod;


%%

void
matrix_error(std::string s)
{
  mparser->error(std::move(s));
}


