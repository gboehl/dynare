/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/KronVector.h,v 1.1.1.1 2004/06/04 13:00:31 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef KRON_VECTOR_H
#define KRON_VECTOR_H

#include "Vector.hh"

class ConstKronVector;

class KronVector : public Vector
{
protected:
  int m{0};
  int n{0};
  int depth{0};
public:
  KronVector() = default;
  KronVector(int mm, int nn, int dp); // new instance
  KronVector(Vector &v, int mm, int nn, int dp); // conversion
  KronVector(KronVector &, int i); // picks i-th subvector
  // We don't want implict conversion from ConstKronVector, since it's expensive
  explicit KronVector(const ConstKronVector &v); // new instance and copy
  KronVector &operator=(const KronVector &v) = default;
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
protected:
  int m;
  int n;
  int depth;
public:
  // Implicit conversion from KronVector is ok, since it's cheap
  ConstKronVector(const KronVector &v);
  ConstKronVector(const ConstKronVector &v);
  ConstKronVector(const Vector &v, int mm, int nn, int dp);
  ConstKronVector(ConstVector v, int mm, int nn, int dp);
  ConstKronVector(const KronVector &v, int i);
  ConstKronVector(const ConstKronVector &v, int i);
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

int power(int m, int depth);

#endif /* KRON_VECTOR */
