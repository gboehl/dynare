// Copyright © 2006, Ondra Kamenik

// $Id: dynare_exception.h 853 2006-08-01 08:42:42Z kamenik $

#ifndef DYNARE_EXCEPTION_H
#define DYNARE_EXCEPTION_H

#include <string>
#include <utility>

class DynareException
{
  std::string mes;
public:
  DynareException(const std::string &m, const std::string &fname, int line, int col)
    : mes{"Parse error at " + fname + ", line " + std::to_string(line) + ", column " +
          std::to_string(col) + ": " + m}
  {
  }
  DynareException(const std::string &fname, int line, const std::string &m)
    : mes{fname + ':' + std::to_string(line) + ": " + m}
  {
  }
  DynareException(const std::string &m, int offset)
    : mes{"Parse error in provided string at offset " + std::to_string(offset) + ": " + m}
  {
  }
  const std::string &
  message() const
  {
    return mes;
  }
};

#endif
