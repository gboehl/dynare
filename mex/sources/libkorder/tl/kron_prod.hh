/*
 * Copyright © 2004 Ondra Kamenik
 * Copyright © 2019-2024 Dynare Team
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

// Kronecker product.

/* Here we define an abstraction for a Kronecker product of a sequence of
   matrices, i.e. A₁⊗…⊗Aₙ. Obviously we do not store the product in memory.
   First we need to represent a dimension of the Kronecker product. Then we
   represent the Kronecker product, simply it is the Kronecker product
   dimension with a vector of references to the matrices A₁,…,Aₙ.

   The main task of this class is to calculate a matrix product B·(A₁⊗A₂⊗…⊗Aₙ),
   which in our application has much more moderate dimensions than A₁⊗A₂⊗…⊗Aₙ.
   We calculate it as B·(A₁⊗I)·…·(I⊗Aᵢ⊗I)·…·(I⊗Aₙ) where dimensions of identity
   matrices differ and are given by the chosen order. One can naturally ask,
   whether there is some optimal order minimizing maximum storage needed for
   intermediate results. The optimal ordering is implemented by class
   KronProdAllOptim.

   For this multiplication, we also need to represent products of type A⊗I,
   I⊗A⊗I and I⊗A. */

#ifndef KRON_PROD_HH
#define KRON_PROD_HH

#include <memory>
#include <utility>
#include <vector>

#include "int_sequence.hh"
#include "permutation.hh"
#include "twod_matrix.hh"

class KronProdAll;
class KronProdAllOptim;
class KronProdIA;
class KronProdIAI;
class KronProdAI;

/* KronProdDimens maintains a dimension of the Kronecker product. So, it
   maintains two sequences, one for rows, and one for columns. */

class KronProdDimens
{
  friend class KronProdAll;
  friend class KronProdAllOptim;
  friend class KronProdIA;
  friend class KronProdIAI;
  friend class KronProdAI;

private:
  IntSequence rows;
  IntSequence cols;

public:
  // Initializes to a given dimension, and all rows and cols are set to zeros
  KronProdDimens(int dim) : rows(dim, 0), cols(dim, 0)
  {
  }
  KronProdDimens(const KronProdDimens& kd) = default;
  KronProdDimens(KronProdDimens&& kd) = default;

  /* Takes dimensions of A₁⊗…⊗Aₙ, and makes dimensions of A₁⊗I, I⊗Aᵢ⊗I or I⊗Aₙ
     for a given i. The dimensions of identity matrices are such that
     A₁⊗…⊗Aₙ=(A₁⊗I)·…·(I⊗Aᵢ⊗I)·…·(I⊗Aₙ). Note that the matrices on the right do
     not commute only because sizes of identity matrices which are then given
     by this ordering. */
  KronProdDimens(const KronProdDimens& kd, int i);

  KronProdDimens& operator=(const KronProdDimens& kd) = default;
  KronProdDimens& operator=(KronProdDimens&& kd) = default;
  [[nodiscard]] bool
  operator==(const KronProdDimens& kd) const
  {
    return rows == kd.rows && cols == kd.cols;
  }

  [[nodiscard]] int
  dimen() const
  {
    return rows.size();
  }
  void
  setRC(int i, int r, int c)
  {
    rows[i] = r;
    cols[i] = c;
  }
  [[nodiscard]] std::pair<int, int>
  getRC(int i) const
  {
    return {rows[i], cols[i]};
  }
  [[nodiscard]] std::pair<int, int>
  getRC() const
  {
    return {rows.mult(), cols.mult()};
  }
  [[nodiscard]] int
  nrows() const
  {
    return rows.mult();
  }
  [[nodiscard]] int
  ncols() const
  {
    return cols.mult();
  }
  [[nodiscard]] int
  nrows(int i) const
  {
    return rows[i];
  }
  [[nodiscard]] int
  ncols(int i) const
  {
    return cols[i];
  }
};

/* Here we define an abstract class for all Kronecker product classes, which
   are KronProdAll (the most general), KronProdIA (for I⊗A), KronProdAI (for
   A⊗I), and KronProdIAI (for I⊗A⊗I). The purpose of the super class is to only
   define some common methods and common member ‘kpd’ for dimensions and
   declare pure virtual mult() which is implemented by the subclasses.

   The class also contains a static method kronMult(), which calculates a
   Kronecker product of two vectors and stores it in the provided vector. It is
   useful at a few points of the library. */

class KronProd
{
protected:
  KronProdDimens kpd;

public:
  KronProd(int dim) : kpd(dim)
  {
  }
  KronProd(KronProdDimens kd) : kpd {std::move(kd)}
  {
  }
  KronProd(const KronProd& kp) = default;
  KronProd(KronProd&& kp) = default;
  virtual ~KronProd() = default;

  [[nodiscard]] int
  dimen() const
  {
    return kpd.dimen();
  }

  virtual void mult(const ConstTwoDMatrix& in, TwoDMatrix& out) const = 0;
  void
  mult(const TwoDMatrix& in, TwoDMatrix& out) const
  {
    mult(ConstTwoDMatrix(in), out);
  }

  void checkDimForMult(const ConstTwoDMatrix& in, const TwoDMatrix& out) const;
  void
  checkDimForMult(const TwoDMatrix& in, const TwoDMatrix& out) const
  {
    checkDimForMult(ConstTwoDMatrix(in), out);
  }

  static void kronMult(const ConstVector& v1, const ConstVector& v2, Vector& res);

  [[nodiscard]] int
  nrows() const
  {
    return kpd.nrows();
  }
  [[nodiscard]] int
  ncols() const
  {
    return kpd.ncols();
  }
  [[nodiscard]] int
  nrows(int i) const
  {
    return kpd.nrows(i);
  }
  [[nodiscard]] int
  ncols(int i) const
  {
    return kpd.ncols(i);
  }
};

/* KronProdAll is the main class of this file. It represents the Kronecker
   product A₁⊗A₂⊗…⊗Aₙ. Besides dimensions, it stores pointers to matrices in
   ‘matlist’ array. If a pointer is null, then the matrix is considered to be
   unit. The array is set by calls to setMat() method (for real matrices) or
   setUnit() method (for unit matrices).

   The object is constructed by a constructor, which allocates the ‘matlist’
   and initializes dimensions to zeros. Then a caller must feed the object with
   matrices by calling setMat() and setUnit() repeatedly for different indices.

   We implement the mult() method of KronProd, and a new method multRows(),
   which creates a vector of kronecker product of all rows of matrices in the
   object. The rows are given by the IntSequence. */

class KronProdAll : public KronProd
{
  friend class KronProdIA;
  friend class KronProdIAI;
  friend class KronProdAI;

protected:
  std::vector<const TwoDMatrix*> matlist;

public:
  KronProdAll(int dim) : KronProd(dim), matlist(dim)
  {
  }
  ~KronProdAll() override = default;
  void setMat(int i, const TwoDMatrix& m);
  void setUnit(int i, int n);
  [[nodiscard]] const TwoDMatrix&
  getMat(int i) const
  {
    return *(matlist[i]);
  }

  void mult(const ConstTwoDMatrix& in, TwoDMatrix& out) const override;
  [[nodiscard]] std::unique_ptr<Vector> multRows(const IntSequence& irows) const;

private:
  [[nodiscard]] bool isUnit() const;
};

/* The class KronProdAllOptim minimizes memory consumption of the product
   B·(A₁⊗A₂⊗…⊗Aₖ). The optimization is done by reordering of the matrices
   A₁,…,Aₖ, in order to minimize a sum of all storages needed for intermediate
   results. The optimal ordering is also nearly optimal with respect to number
   of flops.

   Let (mᵢ,nᵢ) be dimensions of Aᵢ. It is easy to observe, that for the i-th
   step we need storage of r·(n₁·…·nᵢ·mᵢ₊₁·…·mₖ), where r is the number of rows
   of B. To minimize the sum through all i over all permutations of matrices,
   it is equivalent to minimizing the sum
    ₖ
    ∑ (mᵢ₊₁·…·mₖ)/(nᵢ₊₁·…·nₖ)
   ⁱ⁼¹
   The optimal ordering will yield mₖ/nₖ ≤ mₖ₋₁/nₖ₋₁ ≤ … ≤ m₁/n₁

   Now observe, that the number of flops for the i-th step is
   r·(n₁·…·nᵢ·mᵢ·…·mₖ). In order to
   minimize the number of flops, it is equivalent to minimize
    ₖ
    ∑ mᵢ(mᵢ₊₁·…·mₖ)/(nᵢ₊₁·…·nₖ)
   ⁱ⁼¹
   Note that, normally, mᵢ does not change as much as nᵢ₊₁,…,nₖ, so the
   ordering minimizing the memory will be nearly optimal with respect to number
   of flops.

   The class KronProdAllOptim inherits from KronProdAll. A public method
   optimizeOrder() does the reordering. The permutation is stored in ‘oper’. So
   as long as optimizeOrder() is not called, the class is equivalent to
   KronProdAll. */

class KronProdAllOptim : public KronProdAll
{
protected:
  Permutation oper;

public:
  KronProdAllOptim(int dim) : KronProdAll(dim), oper(dim)
  {
  }
  void optimizeOrder();
  [[nodiscard]] const Permutation&
  getPer() const
  {
    return oper;
  }
};

/* This class represents I⊗A. We have only one reference to the matrix, which
   is set by constructor. */

class KronProdIA : public KronProd
{
  friend class KronProdAll;
  const TwoDMatrix& mat;

public:
  KronProdIA(const KronProdAll& kpa) :
      KronProd(KronProdDimens(kpa.kpd, kpa.dimen() - 1)), mat(kpa.getMat(kpa.dimen() - 1))
  {
  }
  void mult(const ConstTwoDMatrix& in, TwoDMatrix& out) const override;
};

/* This class represents A⊗I. We have only one reference to the matrix, which
   is set by constructor. */

class KronProdAI : public KronProd
{
  friend class KronProdIAI;
  friend class KronProdAll;
  const TwoDMatrix& mat;

public:
  KronProdAI(const KronProdAll& kpa) : KronProd(KronProdDimens(kpa.kpd, 0)), mat(kpa.getMat(0))
  {
  }
  KronProdAI(const KronProdIAI& kpiai);

  void mult(const ConstTwoDMatrix& in, TwoDMatrix& out) const override;
};

/* This class represents I⊗A⊗I. We have only one reference to the matrix, which
   is set by constructor. */

class KronProdIAI : public KronProd
{
  friend class KronProdAI;
  friend class KronProdAll;
  const TwoDMatrix& mat;

public:
  KronProdIAI(const KronProdAll& kpa, int i) :
      KronProd(KronProdDimens(kpa.kpd, i)), mat(kpa.getMat(i))
  {
  }
  void mult(const ConstTwoDMatrix& in, TwoDMatrix& out) const override;
};

#endif
