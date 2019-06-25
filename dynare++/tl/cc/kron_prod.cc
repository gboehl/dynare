/*
 * Copyright © 2004 Ondra Kamenik
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

#include "kron_prod.hh"
#include "tl_exception.hh"

#include <tuple>

/* Here we construct Kronecker product dimensions from Kronecker product
   dimensions by picking a given matrix and all other set to identity. The
   constructor takes dimensions of A₁⊗A₂⊗…⊗Aₙ and makes dimensions of A₁⊗I,
   I⊗Aᵢ⊗I or I⊗Aₙ for a given i. The identity matrices must fit into the
   described order. See header file.

   We first decide what is the length of the resulting dimensions. Possible
   length is three for I⊗A⊗I, and two for I⊗A or A⊗I.

   Then we fork according to i. */

KronProdDimens::KronProdDimens(const KronProdDimens &kd, int i)
  : rows((i == 0 || i == kd.dimen()-1) ? 2 : 3),
    cols((i == 0 || i == kd.dimen()-1) ? 2 : 3)
{
  TL_RAISE_IF(i < 0 || i >= kd.dimen(),
              "Wrong index for pickup in KronProdDimens constructor");

  int kdim = kd.dimen();
  if (i == 0)
    {
      // set A⊗I dimensions
      /* The first rows and cols are taken from ‘kd’. The dimension of the
         identity matrix is the number of rows in A₂⊗…⊗Aₙ since the matrix A₁⊗I
         is the first. */
      rows[0] = kd.rows[0];
      rows[1] = kd.rows.mult(1, kdim);
      cols[0] = kd.cols[0];
      cols[1] = rows[1];
    }
  else if (i == kdim-1)
    {
      // set I⊗A dimensions
      /* The second dimension is taken from ‘kd’. The dimension of the identity
         matrix is the number of columns of A₁⊗…⊗Aₙ₋₁, since the matrix I⊗Aₙ is
         the last. */
      rows[0] = kd.cols.mult(0, kdim-1);
      rows[1] = kd.rows[kdim-1];
      cols[0] = rows[0];
      cols[1] = kd.cols[kdim-1];
    }
  else
    {
      // set I⊗A⊗I dimensions
      /* The dimensions of the middle matrix are taken from ‘kd’. The dimension
         of the first identity matrix is the number of columns of A₁⊗…⊗Aᵢ₋₁,
         and the dimension of the last identity matrix is the number of rows of
         Aᵢ₊₁⊗…⊗Aₙ. */
      rows[0] = kd.cols.mult(0, i);
      cols[0] = rows[0];
      rows[1] = kd.rows[i];
      cols[1] = kd.cols[i];
      cols[2] = kd.rows.mult(i+1, kdim);
      rows[2] = cols[2];
    }
}

/* This raises an exception if dimensions are bad for multiplication
   out = in·this. */

void
KronProd::checkDimForMult(const ConstTwoDMatrix &in, const TwoDMatrix &out) const
{
  int my_rows, my_cols;
  std::tie(my_rows, my_cols) = kpd.getRC();
  TL_RAISE_IF(in.nrows() != out.nrows() || in.ncols() != my_rows,
              "Wrong dimensions for KronProd in KronProd::checkDimForMult");
}

/* Here we Kronecker multiply two given vectors v₁ and v₂ and
   store the result in preallocated ‘res’. */

void
KronProd::kronMult(const ConstVector &v1, const ConstVector &v2,
                   Vector &res)
{
  TL_RAISE_IF(res.length() != v1.length()*v2.length(),
              "Wrong vector lengths in KronProd::kronMult");
  res.zeros();
  for (int i = 0; i < v1.length(); i++)
    {
      Vector sub(res, i *v2.length(), v2.length());
      sub.add(v1[i], v2);
    }
}

void
KronProdAll::setMat(int i, const TwoDMatrix &m)
{
  matlist[i] = &m;
  kpd.setRC(i, m.nrows(), m.ncols());
}

void
KronProdAll::setUnit(int i, int n)
{
  matlist[i] = nullptr;
  kpd.setRC(i, n, n);
}

bool
KronProdAll::isUnit() const
{
  int i = 0;
  while (i < dimen() && matlist[i] == nullptr)
    i++;
  return i == dimen();
}

/* Here we compute B·(I⊗A). If m is the dimension of the identity matrix, and B
   is partitioned accordingly, then the result is (B₁·A B₂·A … Bₘ·A).

   In the implementation, ‘outi’ are partitions of ‘out’, ‘ini’ are const
   partitions of ‘in’, and ‘id_cols’ is m. We employ level-2 BLAS. */

void
KronProdIA::mult(const ConstTwoDMatrix &in, TwoDMatrix &out) const
{
  checkDimForMult(in, out);

  int id_cols = kpd.cols[0];
  ConstTwoDMatrix a(mat);

  for (int i = 0; i < id_cols; i++)
    {
      TwoDMatrix outi(out, i * a.ncols(), a.ncols());
      ConstTwoDMatrix ini(in, i * a.nrows(), a.nrows());
      outi.mult(ini, a);
    }
}

/* Here we construct KronProdAI from KronProdIAI. It is clear. */
KronProdAI::KronProdAI(const KronProdIAI &kpiai)
  : KronProd(KronProdDimens(2)), mat(kpiai.mat)
{
  kpd.rows[0] = mat.nrows();
  kpd.cols[0] = mat.ncols();
  kpd.rows[1] = kpiai.kpd.rows[2];
  kpd.cols[1] = kpiai.kpd.cols[2];
}

/* Here we compute B·(A⊗I). Let the dimension of the
   matrix A be m×n, the dimension of I be p, and the number
   of rows of B be q. We use the fact that:

    B·(A⊗I)=reshape(reshape(B, q, mp)·A, q, np).

   This works only for matrix B, whose storage has leading dimension equal to
   number of rows.

   For cases where the leading dimension is not equal to the number of
   rows, we partition the matrix A⊗I into m×n square partitions aᵢⱼI.
   Therefore, we partition B into m partitions (B₁ B₂ … Bₘ). Each partition of B has the same number of
   columns as the identity matrix. If R denotes the resulting matrix,
   then it can be partitioned into n partitions
   (R₁ R₂ … Rₙ). Each partition of R has the same number of
   columns as the identity matrix. Then we have Rᵢ=∑aⱼᵢBⱼ.

   In the implementation, ‘outi’ is Rᵢ, ‘ini’ is Bⱼ, and ‘id_cols’ is the
   dimension of the identity matrix. */

void
KronProdAI::mult(const ConstTwoDMatrix &in, TwoDMatrix &out) const
{
  checkDimForMult(in, out);

  int id_cols = kpd.cols[1];
  ConstTwoDMatrix a(mat);

  if (in.getLD() == in.nrows())
    {
      ConstTwoDMatrix in_resh(in.nrows()*id_cols, a.nrows(), in.getData());
      TwoDMatrix out_resh(in.nrows()*id_cols, a.ncols(), out.getData());
      out_resh.mult(in_resh, a);
    }
  else
    {
      out.zeros();
      for (int i = 0; i < a.ncols(); i++)
        {
          TwoDMatrix outi(out, i *id_cols, id_cols);
          for (int j = 0; j < a.nrows(); j++)
            {
              ConstTwoDMatrix ini(in, j *id_cols, id_cols);
              outi.add(a.get(j, i), ini);
            }
        }
    }
}

/* Here we compute R=B·(I⊗A⊗I). If n is the dimension of the first identity
   matrix, and B is partitioned accordingly, then the result is:
    (B₁·(A⊗I) B₂·(A⊗I) … Bₘ·(A⊗I)).
   Each Bᵢ·(A⊗I) is in fact KronProdAI::mult(). Note that the number of columns
   of partitions of B is the number of rows of A⊗I, and the number of columns
   of partitions of R is the number of columns of A⊗I.

   In the implementation, ‘id_cols’ is n, ‘akronid’ is A⊗I, and ‘in_bl_width’
   and ‘out_bl_width’ are the rows and cols of A⊗I. */

void
KronProdIAI::mult(const ConstTwoDMatrix &in, TwoDMatrix &out) const
{
  checkDimForMult(in, out);

  int id_cols = kpd.cols[0];

  KronProdAI akronid(*this);
  int in_bl_width, out_bl_width;
  std::tie(in_bl_width, out_bl_width) = akronid.kpd.getRC();

  for (int i = 0; i < id_cols; i++)
    {
      TwoDMatrix outi(out, i *out_bl_width, out_bl_width);
      ConstTwoDMatrix ini(in, i *in_bl_width, in_bl_width);
      akronid.mult(ini, outi);
    }
}

/* Here we compute B·(A₁⊗…⊗Aₙ). First we compute B·(A₁⊗I), then this is
   multiplied by all I⊗Aᵢ⊗I, and finally by I⊗Aₙ.

   If the dimension of the Kronecker product is only 1, then we multiply
   two matrices in a straight manner and return.

   The intermediate results are stored on heap pointed by ‘last’. A new
   result is allocated, and then the former storage is deallocated.

   We have to be careful in cases when last or first matrix is unit and
   no calculations are performed in corresponding codes. The codes should
   handle ‘last’ safely also if no calcs are done. */

void
KronProdAll::mult(const ConstTwoDMatrix &in, TwoDMatrix &out) const
{
  // quick copy if product is unit
  if (isUnit())
    {
      out.zeros();
      out.add(1.0, in);
      return;
    }

  // quick zero if one of the matrices is zero
  /* If one of the matrices is exactly zero or the ‘in’ matrix is zero,
     set out to zero and return */
  bool is_zero = false;
  for (int i = 0; i < dimen() && !is_zero; i++)
    is_zero = matlist[i] && matlist[i]->isZero();
  if (is_zero || in.isZero())
    {
      out.zeros();
      return;
    }

  // quick multiplication if dimension is 1
  if (dimen() == 1)
    {
      if (matlist[0]) // always true
        out.mult(in, ConstTwoDMatrix(*(matlist[0])));
      return;
    }

  int c;
  std::unique_ptr<TwoDMatrix> last;

  // perform first multiplication by A₁⊗I
  /* Here we have to construct A₁⊗I, allocate intermediate result ‘last’, and
     perform the multiplication. */
  if (matlist[0])
    {
      KronProdAI akronid(*this);
      c = akronid.kpd.ncols();
      last = std::make_unique<TwoDMatrix>(in.nrows(), c);
      akronid.mult(in, *last);
    }
  else
    last = std::make_unique<TwoDMatrix>(in.nrows(), in.ncols(), Vector{in.getData()});

  // perform intermediate multiplications by I⊗Aᵢ⊗I
  /* Here we go through all I⊗Aᵢ⊗I, construct the product, allocate new storage
     for result ‘newlast’, perform the multiplication and set ‘last’ to
     ‘newlast’. */
  for (int i = 1; i < dimen()-1; i++)
    if (matlist[i])
      {
        KronProdIAI interkron(*this, i);
        c = interkron.kpd.ncols();
        auto newlast = std::make_unique<TwoDMatrix>(in.nrows(), c);
        interkron.mult(*last, *newlast);
        last = std::move(newlast);
      }

  // perform last multiplication by I⊗Aₙ
  if (matlist[dimen()-1])
    {
      KronProdIA idkrona(*this);
      idkrona.mult(*last, out);
    }
  else
    out = *last;
}

/* This calculates a Kornecker product of rows of matrices, the row
   indices are given by the integer sequence. The result is allocated and
   returned. */

std::unique_ptr<Vector>
KronProdAll::multRows(const IntSequence &irows) const
{
  TL_RAISE_IF(irows.size() != dimen(),
              "Wrong length of row indices in KronProdAll::multRows");

  std::unique_ptr<Vector> last;
  std::unique_ptr<ConstVector> row;
  std::vector<std::unique_ptr<Vector>> to_delete;
  for (int i = 0; i < dimen(); i++)
    {
      int j = dimen()-1-i;

      // set ‘row’ to the number of rows of j-th matrix
      /* If the j-th matrix is a real matrix, then the row is constructed from
         the matrix. It the matrix is unit, we construct a new vector, fill it
         with zeros, than set the unit to appropriate place, and make the ‘row’
         as ConstVector of this vector, which is sheduled for deallocation. */
      if (matlist[j])
        row = std::make_unique<ConstVector>(matlist[j]->getRow(irows[j]));
      else
        {
          auto aux = std::make_unique<Vector>(ncols(j));
          aux->zeros();
          (*aux)[irows[j]] = 1.0;
          row = std::make_unique<ConstVector>(*aux);
          to_delete.emplace_back(std::move(aux));
        }

      // set ‘last’ to product of ‘row’ and ‘last’
      /* If the ‘last’ exists, we allocate new storage and Kronecker multiply.
         If the ‘last’ does not exist, then we only make ‘last’ equal to
         ‘row’. */
      if (last)
        {
          auto newlast = std::make_unique<Vector>(last->length()*row->length());
          kronMult(*row, ConstVector(*last), *newlast);
          last = std::move(newlast);
        }
      else
        last = std::make_unique<Vector>(*row);
    }

  return last;
}

/* This permutes the matrices so that the new ordering would minimize memory
   consumption. As shown in KronProdAllOptim class declaration, we want:

     mₖ/nₖ ≤ mₖ₋₁/nₖ₋₁ ≤ … ≤ m₁/n₁

   where mᵢ×nᵢ is the dimension of Aᵢ. So we implement the bubble sort. */

void
KronProdAllOptim::optimizeOrder()
{
  for (int i = 0; i < dimen(); i++)
    {
      int swaps = 0;
      for (int j = 0; j < dimen()-1; j++)
        {
          if (static_cast<double>(kpd.rows[j])/kpd.cols[j]
              < static_cast<double>(kpd.rows[j+1])/kpd.cols[j+1])
            {
              // swap dimensions and matrices at j and j+1
              int s = kpd.rows[j+1];
              kpd.rows[j+1] = kpd.rows[j];
              kpd.rows[j] = s;
              s = kpd.cols[j+1];
              kpd.cols[j+1] = kpd.cols[j];
              kpd.cols[j] = s;
              const TwoDMatrix *m = matlist[j+1];
              matlist[j+1] = matlist[j];
              matlist[j] = m;

              // project the swap to the permutation ‘oper’
              s = oper.getMap()[j+1];
              oper.getMap()[j+1] = oper.getMap()[j];
              oper.getMap()[j] = s;
              swaps++;
            }
        }
      if (swaps == 0)
        return;
    }
}
