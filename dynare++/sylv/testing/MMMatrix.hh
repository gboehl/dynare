/* $Header: /var/lib/cvs/dynare_cpp/sylv/testing/MMMatrix.h,v 1.1.1.1 2004/06/04 13:01:13 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef MM_MATRIX_H
#define MM_MATRIX_H

#include "GeneralMatrix.hh"
#include "SylvMemory.hh"

#include <string>
#include <utility>
#include <memory>

class MMException : public MallocAllocator
{
  std::string message;
public:
  MMException(std::string mes) : message(std::move(mes))
  {
  }
  std::string
  getMessage() const
  {
    return message;
  }
};

class MMMatrixIn : public MallocAllocator
{
  std::shared_ptr<double> data;
  int rows;
  int cols;
public:
  MMMatrixIn(const std::string &fname);
  ~MMMatrixIn() = default;
  Vector
  getData() const
  {
    return Vector{data, size()};
  }
  int
  size() const
  {
    return rows*cols;
  }
  int
  row() const
  {
    return rows;
  }
  int
  col() const
  {
    return cols;
  }
};

class MMMatrixOut : public MallocAllocator
{
public:
  static void write(const std::string &fname, const GeneralMatrix &m);
};

#endif /* MM_MATRIX_H */
