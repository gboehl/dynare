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
#include "formula_tab.hh"

#define YY_USER_ACTION SET_LLOC(fmla_);
%}

%option nounput
%option noyy_top_state
%option stack
%option yylineno
%option prefix="fmla_"
%option never-interactive
%x CMT

%%

 /* comments */
<*>"/*"              {yy_push_state(CMT);}
<CMT>[^*\n]*
<CMT>"*"+[^*/\n]*
<CMT>"*"+"/"         {yy_pop_state();}
<CMT>[\n]
"//".*\n

 /* initial spaces or tabs are ignored */

[ \t\r\n]
[+]                  {return YPLUS;}
[-]                  {return YMINUS;}
[*]                  {return YTIMES;}
[/]                  {return YDIVIDE;}
[\^]                 {return YPOWER;}
exp                  {return YEXP;}
log                  {return YLOG;}
sin                  {return YSIN;}
cos                  {return YCOS;}
tan                  {return YTAN;}
sqrt                 {return YSQRT;}
erf                  {return YERF;}
erfc                 {return YERFC;}
diff                 {return YDIFF;}

 /* names: parameters, variables (lagged/leaded) */
[A-Za-z_][A-Za-z0-9_]*([\(\{][+-]?[0-9]+[\)\}])? {
	fmla_lval.string=fmla_text;
	return NAME;
}

 /* floating point numbers */
(([0-9]*\.?[0-9]+)|([0-9]+\.))([edED][-+]?[0-9]+)? {
	fmla_lval.string=fmla_text;
	return DNUMBER;
}

=                    {return EQUAL_SIGN;}

.                    {return fmla_text[0];}

%%

int
fmla_wrap()
{
  return 1;
}

void
fmla__destroy_buffer(void* p)
{
  fmla__delete_buffer(static_cast<YY_BUFFER_STATE>(p));
}
