// Copyright 2005, Ondra Kamenik

// Exception

/* This is a simple code defining an exception and two convenience macros. */

#include <string>
#include <iostream>

#ifndef KORD_EXCEPTION_H
#define KORD_EXCEPTION_H

#define KORD_RAISE(mes)                         \
  throw KordException(__FILE__, __LINE__, mes);

#define KORD_RAISE_IF(expr, mes)                                \
  if (expr) throw KordException(__FILE__, __LINE__, mes);

#define KORD_RAISE_X(mes, c)                            \
  throw KordException(__FILE__, __LINE__, mes, c);

#define KORD_RAISE_IF_X(expr, mes, c)                           \
  if (expr) throw KordException(__FILE__, __LINE__, mes, c);

class KordException
{
protected:
  std::string fname;
  int lnum;
  std::string message;
  int cd;
public:
  KordException(std::string f, int l, std::string mes, int c = 255)
    : fname{std::move(f)}, lnum{l}, message{std::move(mes)}, cd{c}
  {
  }
  virtual ~KordException() = default;
  virtual void
  print() const
  {
    std::cout << "At " << fname << ':' << lnum << ":(" << cd << "):" << message << '\n';
  }
  virtual int
  code() const
  {
    return cd;
  }
  const std::string &
  get_message() const
  {
    return message;
  }
};

// |KordException| error code definitions
constexpr int KORD_FP_NOT_CONV = 254;
constexpr int KORD_FP_NOT_FINITE = 253;
constexpr int KORD_MD_NOT_STABLE = 252;

#endif
