/* $Header: /var/lib/cvs/dynare_cpp/sylv/testing/MMMatrix.h,v 1.1.1.1 2004/06/04 13:01:13 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef MM_MATRIX_H
#define MM_MATRIX_H

#include "GeneralMatrix.hh"
#include "SylvMemory.hh"

#include <string>
#include <utility>
#include <memory>

using namespace std;

class MMException : public MallocAllocator
{
  string message;
public:
  MMException(string mes) : message(std::move(mes))
  {
  }
  MMException(const char *mes) : message(mes)
  {
  }
  const char *
  getMessage() const
  {
    return message.data();
  }
};

class MMMatrixIn : public MallocAllocator
{
  std::shared_ptr<const double> data;
  int rows;
  int cols;
public:
  MMMatrixIn(const char *fname);
  ~MMMatrixIn() = default;
  ConstVector
  getData() const
  {
    return ConstVector{data, size()};
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
  static void write(const char *fname, int rows, int cols, const double *data);
  static void write(const char *fname, const GeneralMatrix &m);
};

#endif /* MM_MATRIX_H */
