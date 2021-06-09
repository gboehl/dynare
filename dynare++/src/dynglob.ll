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
#include "parser/cc/location.hh"
#include "dynglob_tab.hh"

#define YY_USER_ACTION SET_LLOC(dynglob_);
%}

%option nounput
%option noyy_top_state
%option stack
%option yylineno
%option prefix="dynglob_"
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

 /* initial spaces or tabs are ignored */

[ \t\r\n\0]
var                {return VAR;}
varexo             {return VAREXO;}
parameters         {return PARAMETERS;}
model              {return MODEL;}
end                {return END;}
initval            {return INITVAL;}
order              {return ORDER;}
vcov               {return VCOV;}
planner_objective  {return PLANNEROBJECTIVE;}
planner_discount   {return PLANNERDISCOUNT;}

 /* names */
[A-Za-z_][A-Za-z0-9_]* {
	dynglob_lval.string = dynglob_text;
	return NAME;
}

;                  {return SEMICOLON;}
,                  {return COMMA;}
=                  {return EQUAL_SIGN;}
\[                 {return LEFT_BRACKET;}
\]                 {return RIGHT_BRACKET;}
. {
	dynglob_lval.character = dynglob_text[0];
	return CHARACTER;
}

%%

int
dynglob_wrap()
{
  return 1;
}

void
dynglob__destroy_buffer(void* p)
{
  dynglob__delete_buffer(static_cast<YY_BUFFER_STATE>(p));
}
