// Copyright (C) 2005, Ondra Kamenik

// $Id: pascal_triangle.h 762 2006-05-22 13:00:07Z kamenik $

#ifndef PASCAL_TRIANGLE_H
#define PASCAL_TRIANGLE_H

#include <vector>

class PascalRow : public std::vector<int>
{
  int k{1};
public:
  PascalRow() : std::vector<int>{}
  {
    push_back(2);
  }
  void setFromPrevious(const PascalRow &prev);
  void prolong(const PascalRow &prev);
  void prolongFirst(int n);
  void print() const;
};

namespace PascalTriangle
{
  int noverk(int n, int k);
  void print();
};

#endif
