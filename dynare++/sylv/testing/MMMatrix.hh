/* $Header: /var/lib/cvs/dynare_cpp/sylv/testing/MMMatrix.h,v 1.1.1.1 2004/06/04 13:01:13 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef MM_MATRIX_H
#define MM_MATRIX_H

#include "GeneralMatrix.hh"

#include <string>
#include <utility>
#include <vector>

class MMException
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

class MMMatrixIn
{
  std::vector<double> data;
  int rows;
  int cols;
public:
  MMMatrixIn(const std::string &fname);
  ~MMMatrixIn() = default;
  Vector
  getData()
  {
    return Vector{data.data(), size()};
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

class MMMatrixOut
{
public:
  static void write(const std::string &fname, const GeneralMatrix &m);
};

#endif /* MM_MATRIX_H */
