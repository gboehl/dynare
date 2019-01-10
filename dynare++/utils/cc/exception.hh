// Copyright (C) 2005, Ondra Kamenik

// $Id: exception.h 1367 2007-07-11 14:21:57Z kamenik $

#ifndef OGU_EXCEPTION_H
#define OGU_EXCEPTION_H

#include <string>
#include <iostream>
#include <utility>

namespace ogu
{
  /** A primitive exception. */
  class Exception
  {
  protected:
    const std::string file;
    const int line;
    const std::string mes;
  public:
    Exception(std::string file_arg, int line_arg, std::string mes_arg)
      : file{std::move(file_arg)},
        line{line_arg},
        mes{std::move(mes_arg)}
    {
    }
    virtual ~Exception() = default;

    void
    print(std::ostream &out) const
    {
      out << file << ':' << line << ": " << mes << std::endl;
    }

    void
    print() const
    {
      print(std::cout);
    }

    std::string
    message() const
    {
      return mes;
    }
  };
};

#endif
