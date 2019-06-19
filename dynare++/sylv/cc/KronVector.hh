/*
 * Copyright © 2004-2011 Ondra Kamenik
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

#ifndef KRON_VECTOR_H
#define KRON_VECTOR_H

#include "Vector.hh"

class ConstKronVector;

class KronVector : public Vector
{
  friend class ConstKronVector;
protected:
  int m{0};
  int n{0};
  int depth{0};
public:
  KronVector() = default;
  KronVector(const KronVector &v) = default;
  KronVector(KronVector &&v) = default;
  KronVector(int mm, int nn, int dp); // new instance
  KronVector(Vector &v, int mm, int nn, int dp); // conversion
  KronVector(KronVector &, int i); // picks i-th subvector
  // We don't want implict conversion from ConstKronVector, since it's expensive
  explicit KronVector(const ConstKronVector &v); // new instance and copy
  KronVector &operator=(const KronVector &v) = default;
  KronVector &operator=(KronVector &&v) = default;
  KronVector &operator=(const ConstKronVector &v);
  KronVector &operator=(const Vector &v);
  int
  getM() const
  {
    return m;
  }
  int
  getN() const
  {
    return n;
  }
  int
  getDepth() const
  {
    return depth;
  }
};

class ConstKronVector : public ConstVector
{
  friend class KronVector;
protected:
  int m;
  int n;
  int depth;
public:
  // Implicit conversion from KronVector is ok, since it's cheap
  ConstKronVector(const KronVector &v);
  ConstKronVector(const ConstKronVector &v) = default;
  ConstKronVector(ConstKronVector &&v) = default;
  ConstKronVector(const Vector &v, int mm, int nn, int dp);
  ConstKronVector(ConstVector v, int mm, int nn, int dp);
  ConstKronVector(const KronVector &v, int i);
  ConstKronVector(const ConstKronVector &v, int i);
  ConstKronVector &operator=(const ConstKronVector &v) = delete;
  ConstKronVector &operator=(ConstKronVector &&v) = delete;
  int
  getM() const
  {
    return m;
  }
  int
  getN() const
  {
    return n;
  }
  int
  getDepth() const
  {
    return depth;
  }
};

#endif /* KRON_VECTOR */
