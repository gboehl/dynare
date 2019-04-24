// -*- C++ -*-
/* Copyright © 2006-2011, Ondra Kamenik */

%code requires
{
#include "location.hh"
#define FMLA_LTYPE ogp::location_type
}

%code
{
#include "formula_parser.hh" 
#include <string>

void fmla_error(std::string);
int fmla_lex();
extern ogp::FormulaParser* fparser;
}

%union
{
  char* string;
  double dvalue;
  int integer;
}

%token EQUAL_SIGN
%left YPLUS YMINUS
%left YTIMES YDIVIDE
%precedence YUMINUS YUPLUS
%right YPOWER
%token YEXP YLOG YSIN YCOS YTAN YSQRT YERF YERFC YDIFF
%token <string> DNUMBER NAME
%type <integer> expression

%define api.prefix {fmla_}

%locations
%defines
%define parse.error verbose

%%
 root : equation_list
      | expression
                                  {fparser->add_formula($1);}
      ; 

 equation_list : equation_list equation | equation ;

 equation : expression EQUAL_SIGN expression ';' 
                                  {fparser->add_formula(fparser->add_binary(ogp::code_t::MINUS,$1,$3));}
      | expression ';'
                                  {fparser->add_formula($1);}
      ;

  expression : '(' expression ')' { $$ = $2;}
      | expression YPLUS expression {$$=fparser->add_binary(ogp::code_t::PLUS,$1,$3);}
      | expression YMINUS expression {$$=fparser->add_binary(ogp::code_t::MINUS,$1,$3);}
      | expression YTIMES expression {$$=fparser->add_binary(ogp::code_t::TIMES,$1,$3);}
      | expression YDIVIDE expression {$$=fparser->add_binary(ogp::code_t::DIVIDE,$1,$3);}
      | expression YPOWER expression {$$=fparser->add_binary(ogp::code_t::POWER,$1,$3);}
      | YMINUS expression %prec YUMINUS {$$=fparser->add_unary(ogp::code_t::UMINUS,$2);}
      | YPLUS expression %prec YUPLUS {$$ = $2;}
      | YSIN '(' expression ')' {$$=fparser->add_unary(ogp::code_t::SIN,$3);}
      | YCOS '(' expression ')' {$$=fparser->add_unary(ogp::code_t::COS,$3);}
      | YTAN '(' expression ')' {$$=fparser->add_unary(ogp::code_t::TAN,$3);}
      | YEXP '(' expression ')' {$$=fparser->add_unary(ogp::code_t::EXP,$3);}
      | YLOG '(' expression ')' {$$=fparser->add_unary(ogp::code_t::LOG,$3);}
      | YSQRT '(' expression ')' {$$=fparser->add_unary(ogp::code_t::SQRT,$3);}
      | YERF '(' expression ')' {$$=fparser->add_unary(ogp::code_t::ERF,$3);}
      | YERFC '(' expression ')' {$$=fparser->add_unary(ogp::code_t::ERFC,$3);}
      | YDIFF '(' expression ',' NAME ')' {$$=fparser->add_derivative($3, fparser->add_nulary($5));}
      | NAME {$$=fparser->add_nulary($1);}
      | DNUMBER {$$=fparser->add_nulary($1);}
      ;

%%

void
fmla_error(std::string s)
{
  fparser->error(std::move(s));
}
