/*
 * Copyright © 2005 Ondra Kamenik
 * Copyright © 2019 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

// Product quadrature.

/* This file defines a product multidimensional quadrature. If Qₖ$ denotes the
   one dimensional quadrature, then the product quadrature Q of k level and
   dimension d takes the form:

          nₖ      nₖ  
    Qf =  ∑   …   ∑   w_i₁·…·w_{i_d} f(x_i₁,…,x_{i_d})
         i₁=1   i_d=1

   which can be written in terms of the one dimensional quadrature Qₖ as

    Qf=(Qₖ⊗…⊗Qₖ)f

   Here we define the product quadrature iterator prodpit and plug it into
   QuadratureImpl to obtains ProductQuadrature. */

#ifndef PRODUCT_H
#define PRODUCT_H

#include "int_sequence.hh"
#include "vector_function.hh"
#include "quadrature.hh"

/* This defines a product point iterator. We have to maintain the following: a
   pointer to product quadrature in order to know the dimension and the
   underlying one dimensional quadrature, then level, number of points in the
   level, integer sequence of indices, signal, the coordinates of the point and
   the weight.

   The point indices, signal, and point coordinates are implmented as pointers
   in order to allow for empty constructor.

   The constructor prodpit(const ProductQuadrature& q, int j0, int l)
   constructs an iterator pointing to (j0,0,…,0), which is used by begin()
   dictated by QuadratureImpl. */

class ProductQuadrature;

class prodpit
{
protected:
  const ProductQuadrature &prodq;
  int level{0};
  int npoints{0};
  IntSequence jseq;
  bool end_flag{true};
  ParameterSignal sig;
  Vector p;
  double w;
public:
  prodpit() = default;
  prodpit(const ProductQuadrature &q, int j0, int l);
  prodpit(const prodpit &ppit) = default;
  ~prodpit() = default;
  bool operator==(const prodpit &ppit) const;
  bool
  operator!=(const prodpit &ppit) const
  {
    return !operator==(ppit);
  }
  prodpit &operator=(const prodpit &spit) = delete;
  prodpit &operator++();
  const ParameterSignal &
  signal() const
  {
    return sig;
  }
  const Vector &
  point() const
  {
    return p;
  }
  double
  weight() const
  {
    return w;
  }
  void print() const;
protected:
  void setPointAndWeight();
};

/* The product quadrature is just QuadratureImpl with the product iterator
   plugged in. The object is constructed by just giving the underlying one
   dimensional quadrature, and the dimension. The only extra method is
   designLevelForEvals() which for the given maximum number of evaluations (and
   dimension and underlying quadrature from the object) returns a maximum level
   yeilding number of evaluations less than the given number. */

class ProductQuadrature : public QuadratureImpl<prodpit>
{
  friend class prodpit;
  const OneDQuadrature &uquad;
public:
  ProductQuadrature(int d, const OneDQuadrature &uq);
  ~ProductQuadrature() override = default;
  int
  numEvals(int l) const override
  {
    int res = 1;
    for (int i = 0; i < dimen(); i++)
      res *= uquad.numPoints(l);
    return res;
  }
  void designLevelForEvals(int max_eval, int &lev, int &evals) const;
protected:
  prodpit begin(int ti, int tn, int level) const override;
};

#endif
