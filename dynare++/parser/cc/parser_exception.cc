// Copyright © 2006, Ondra Kamenik

// $Id: parser_exception.cpp 2269 2008-11-23 14:33:22Z michel $

#include "parser_exception.hh"

using namespace ogp;

ParserException::ParserException(string m, int offset)
  : mes(std::move(m)), off(offset),
    aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
}

ParserException::ParserException(string m, const char *dum, int i1)
  : mes(std::move(m)), off(0),
    aux_i1(i1), aux_i2(-1), aux_i3(-1)
{
}

ParserException::ParserException(string m, const char *dum, int i1, int i2)
  : mes(std::move(m)), off(0),
    aux_i1(i1), aux_i2(i2), aux_i3(-1)
{
}

ParserException::ParserException(string m, const char *dum, int i1, int i2, int i3)
  : mes(std::move(m)), off(0),
    aux_i1(i1), aux_i2(i2), aux_i3(i3)
{
}

ParserException::ParserException(const ParserException &m, int plus_offset)
  : aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
  copy(m);
  off += plus_offset;
}

ParserException::ParserException(const ParserException &m, const char *dum, int i)
  : aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
  copy(m);
  aux_i3 = m.aux_i2;
  aux_i2 = m.aux_i1;
  aux_i1 = i;
}

ParserException::ParserException(const ParserException &m, const char *dum, int i1, int i2)
  : aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
  copy(m);
  aux_i3 = m.aux_i1;
  aux_i2 = i2;
  aux_i1 = i1;
}

ParserException::ParserException(const ParserException &m, const char *dum, int i1, int i2, int i3)
  : aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
  copy(m);
  aux_i3 = i3;
  aux_i2 = i2;
  aux_i1 = i1;
}

void
ParserException::copy(const ParserException &e)
{
  mes = e.mes;
  off = e.off;
  aux_i1 = e.aux_i1;
  aux_i2 = e.aux_i2;
  aux_i3 = e.aux_i3;
}
