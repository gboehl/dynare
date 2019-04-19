#include "csv_parser.hh"
#include "parser_exception.hh"
#include "location.hh"
#include "csv_tab.hh"

#include <memory>
#include <algorithm>

using namespace ogp;

/** A global symbol for passing info to the CSVParser from
 * csv_parse(). */
CSVParser *csv_parser;

/** The declaration of functions defined in csv_ll.cc and
 * csv_tab.cc generated from csv.lex and csv.y. */
void *csv__scan_buffer(char *, unsigned int);
void csv__destroy_buffer(void *);
int csv_parse();

extern ogp::location_type csv_lloc;

void
CSVParser::csv_error(const char *mes)
{
  throw ParserException(mes, csv_lloc.off);
}

void
CSVParser::csv_parse(int length, const char *str)
{
  // allocate temporary buffer and parse
  auto buffer = std::make_unique<char[]>(length+2);
  std::copy_n(str, length, buffer.get());
  buffer[length] = '\0';
  buffer[length+1] = '\0';
  csv_lloc.off = 0;
  csv_lloc.ll = 0;
  parsed_string = buffer.get();
  void *p = csv__scan_buffer(buffer.get(), static_cast<unsigned int>(length)+2);
  csv_parser = this;
  ::csv_parse();
  csv__destroy_buffer(p);
  parsed_string = nullptr;
}
