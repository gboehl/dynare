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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

// Smolyak quadrature.

/* This file defines Smolyak (sparse grid) multidimensional quadrature for
   non-nested underlying one dimensional quadrature. Let Q¹ₗ denote the one
   dimensional quadrature of l level. Let nₗ denote a number of points in the l
   level. Than the Smolyak quadrature can be defined as

                                      ⎛ d−1 ⎞
    Qᵈf =     ∑      (−1)^{l+d−|k|−1} ⎢     ⎥ (Q¹_k₁⊗…⊗Q¹_{k_d})f
         l≤|k|≤l+d−1                  ⎝|k|−l⎠  

   where d is the dimension, k is d-dimensional sequence of integers, and |k|
   denotes the sum of the sequence.

   Here we define smolpit as Smolyak iterator and SmolyakQuadrature. */

#ifndef SMOLYAK_H
#define SMOLYAK_H

#include "int_sequence.hh"
#include "tl_static.hh"
#include "vector_function.hh"
#include "quadrature.hh"
#include "pascal_triangle.hh"

/* Here we define the Smolyak point iterator. The Smolyak formula can be broken
   to a sum of product quadratures with various combinations of levels. The
   iterator follows this pattern. It maintains an index to a summand and then a
   point coordinates within the summand (product quadrature). The array of
   summands to which the ‘isummand’ points is maintained by the
   SmolyakQuadrature class to which the object knows the pointer ‘smolq’.

   We provide a constructor which points to the beginning of the given summand.
   This constructor is used in SmolyakQuadrature::begin() method which
   approximately divideds all the iterators to subsets of equal size. */

class SmolyakQuadrature;

class smolpit
{
protected:
  const SmolyakQuadrature &smolq;
  unsigned int isummand{0};
  IntSequence jseq;
  ParameterSignal sig;
  Vector p;
  double w;
public:
  smolpit(const SmolyakQuadrature &q, unsigned int isum);
  smolpit(const smolpit &spit) = default;
  ~smolpit() = default;
  bool operator==(const smolpit &spit) const;
  bool
  operator!=(const smolpit &spit) const
  {
    return !operator==(spit);
  }
  smolpit &operator=(const smolpit &spit) = delete;
  smolpit &operator++();
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

/* Here we define the class SmolyakQuadrature. It maintains an array of
   summands of the Smolyak quadrature formula:

                                 ⎛ d−1 ⎞
         ∑      (−1)^{l+d−|k|−1} ⎢     ⎥ (Q¹_k₁⊗…⊗Q¹_{k_d})f
    l≤|k|≤l+d−1                  ⎝|k|−l⎠  

   Each summand is fully specified by sequence k. The summands are here
   represented (besides k) also by sequence of number of points in each level
   selected by k, and also by a cummulative number of evaluations. The latter
   two are added only for conveniency.

   The summands in the code are given by ‘levels’, which is a vector of
   k sequences, further by ‘levpoints’ which is a vector of sequences
   of nuber of points in each level, and by ‘cumevals’ which is the
   cumulative number of points, this is:

       d
    ∑  ∏  n_kᵢ
    ᵏ ⁱ⁼¹

   where the sum is done through all k before the current.

   The ‘levels’ and ‘levpoints’ vectors are used by smolpit. */

class SmolyakQuadrature : public QuadratureImpl<smolpit>
{
  friend class smolpit;
  int level;
  const OneDQuadrature &uquad;
  std::vector<IntSequence> levels;
  std::vector<IntSequence> levpoints;
  std::vector<int> cumevals;
public:
  SmolyakQuadrature(int d, int l, const OneDQuadrature &uq);
  ~SmolyakQuadrature() override = default;
  int numEvals(int level) const override;
  void designLevelForEvals(int max_eval, int &lev, int &evals) const;
protected:
  smolpit begin(int ti, int tn, int level) const override;
  unsigned int
  numSummands() const
  {
    return levels.size();
  }
private:
  int calcNumEvaluations(int level) const;
};

#endif
