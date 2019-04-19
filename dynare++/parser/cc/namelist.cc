// Copyright © 2006, Ondra Kamenik

// $Id: namelist.cpp 42 2007-01-22 21:53:24Z ondra $

#include "namelist.hh"

#include <memory>
#include <algorithm>

using namespace ogp;

/** A global symbol for passing info to NameListParser from its
 * parser. */
NameListParser *name_list_parser;

void *namelist__scan_buffer(char *, unsigned int);
void namelist__destroy_buffer(void *);
void namelist_parse();

void
NameListParser::namelist_parse(int length, const char *stream)
{
  auto buffer = std::make_unique<char[]>(length+2);
  std::copy_n(str, length, buffer.get());
  buffer[length] = '\0';
  buffer[length+1] = '\0';
  void *p = namelist__scan_buffer(buffer.get(), static_cast<unsigned int>(length)+2);
  name_list_parser = this;
  ::namelist_parse();
  namelist__destroy_buffer(p);
}
