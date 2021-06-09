// -*- C++ -*-
/*
 * Copyright © 2006-2011 Ondra Kamenik
 * Copyright © 2019 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

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


