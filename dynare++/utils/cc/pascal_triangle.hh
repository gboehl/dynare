// Copyright (C) 2005, Ondra Kamenik

// $Id: pascal_triangle.h 762 2006-05-22 13:00:07Z kamenik $

#ifndef PASCAL_TRIANGLE_H
#define PASCAL_TRIANGLE_H

#include <vector>

namespace ogu
{

  using std::vector;

  class PascalRow : public vector<int>
  {
    int k{1};
  public:
    PascalRow()
      : vector<int>() 
    {
      push_back(2);
    }
    void setFromPrevious(const PascalRow &prev);
    void prolong(const PascalRow &prev);
    void prolongFirst(int n);
    void print() const;
  };

  class PascalTriangle
  {
    vector<PascalRow> tr;
  public:
    PascalTriangle()
    {
      tr.emplace_back();
    }
    PascalTriangle(const PascalTriangle &triang)
       
    = default;
    PascalTriangle &
    operator=(const PascalTriangle &triang)
    = default;
    int noverk(int n, int k);
    void print() const;
  protected:
    void ensure(int n, int k);
    int max_n() const;
    int max_k() const;
  };
};

extern ogu::PascalTriangle ptriang;

#endif
