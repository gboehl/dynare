/* $Id: factory.h 148 2005-04-19 15:12:26Z kamenik $ */
/* Copyright 2004, Ondra Kamenik */

#ifndef FACTORY_H
#define FACTORY_H

#include <random>
#include <memory>

#include "symmetry.hh"
#include "int_sequence.hh"
#include "twod_matrix.hh"
#include "equivalence.hh"
#include "rfs_tensor.hh"
#include "t_container.hh"

class Factory
{
  std::mt19937 mtgen;
  std::uniform_real_distribution<> dis;

  void init(const Symmetry &s, const IntSequence &nvs);
  void init(int dim, int nv);
  void fillMatrix(TwoDMatrix &m);
public:
  double get();
  // this can be used with UGSTensor, FGSTensor
  template <class _Ttype>
  std::unique_ptr<_Ttype>
  make(int r, const Symmetry &s, const IntSequence &nvs)
  {
    auto res = std::make_unique<_Ttype>(r, TensorDimens(s, nvs));
    init(s, nvs);
    fillMatrix(*res);
    return res;
  }

  // this can be used with FFSTensor, UFSTensor, FRTensor, URTensor
  template <class _Ttype>
  std::unique_ptr<_Ttype>
  make(int r, int nv, int dim)
  {
    auto res = std::make_unique<_Ttype>(r, nv, dim);
    init(dim, nv);
    fillMatrix(*res);
    return res;
  }

  template <class _Ttype, class _Ctype>
  _Ctype
  makeCont(int r, const IntSequence &nvs, int maxdim)
  {
    int symnum = nvs.size();
    _Ctype res(symnum);
    for (int dim = 1; dim <= maxdim; dim++)
      if (symnum == 1)
        // full symmetry
        res.insert(make<_Ttype>(r, Symmetry{dim}, nvs));
      else
        // general symmetry
        for (int i = 0; i <= dim; i++)
          res.insert(make<_Ttype>(r, Symmetry{i, dim-i}, nvs));
    return res;
  }

  template <class _Ttype, class _Ptype>
  _Ptype
  makePoly(int r, int nv, int maxdim)
  {
    _Ptype p(r, nv);
    for (int d = 1; d <= maxdim; d++)
      p.insert(make<_Ttype>(r, nv, d));
    return p;
  }

  Vector makeVector(int n);
};

#endif
