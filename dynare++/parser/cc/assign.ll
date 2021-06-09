/* -*- C++ -*- */
/*
 * Copyright © 2004-2011 Ondra Kamenik
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
%{
#include "location.hh"
#include "assign_tab.hh"

#define YY_USER_ACTION SET_LLOC(asgn_);
%}

%option nounput
%option noyy_top_state
%option stack
%option yylineno
%option prefix="asgn_"
%option never-interactive
%x CMT

%%

 /* comments */
<*>"/*"            {yy_push_state(CMT);}
<CMT>[^*\n]*
<CMT>"*"+[^*/\n]*
<CMT>"*"+"/"       {yy_pop_state();}
<CMT>[\n]
"//".*\n

 /* spaces */
[ \t\r\n]          {return BLANK;}

 /* names */
[A-Za-z_][A-Za-z0-9_]* {
	asgn_lval.string = asgn_text;
	return NAME;
}

;                  {return SEMICOLON;}
=                  {return EQUAL_SIGN;}
. {
	asgn_lval.character = asgn_text[0];
	return CHARACTER;
}

%%

int
asgn_wrap()
{
  return 1;
}

void
asgn__destroy_buffer(void* p)
{
  asgn__delete_buffer(static_cast<YY_BUFFER_STATE>(p));
}
