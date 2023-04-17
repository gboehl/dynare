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
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "sparse_tensor.hh"
#include "fs_tensor.hh"
#include "tl_exception.hh"

#include <iostream>
#include <iomanip>
#include <cmath>

/* This is straightforward. Before we insert anything, we do a few
   checks. Then we reset ‘first_nz_row’ and ‘last_nz_row’ if necessary. */

void
SparseTensor::insert(IntSequence key, int r, double c)
{
  TL_RAISE_IF(r < 0 || r >= nr,
              "Row number out of dimension of tensor in SparseTensor::insert");
  TL_RAISE_IF(key.size() != dimen(),
              "Wrong length of key in SparseTensor::insert");
  TL_RAISE_IF(!std::isfinite(c),
              "Insertion of non-finite value in SparseTensor::insert");

  auto first_pos = m.lower_bound(key);

  // check that pair ‘key’ and ‘r’ is unique
  auto last_pos = m.upper_bound(key);
  for (auto it = first_pos; it != last_pos; ++it)
    TL_RAISE_IF(it->second.first == r, "Duplicate <key, r> insertion in SparseTensor::insert");

  m.emplace_hint(first_pos, std::move(key), std::make_pair(r, c));
  if (first_nz_row > r)
    first_nz_row = r;
  if (last_nz_row < r)
    last_nz_row = r;
}

/* This returns true if all items are finite (not NaN nor ∞). */

bool
SparseTensor::isFinite() const
{
  bool res = true;
  auto run = m.begin();
  while (res && run != m.end())
    {
      if (!std::isfinite(run->second.second))
        res = false;
      ++run;
    }
  return res;
}

/* This returns a ratio of a number of non-zero columns in folded
   tensor to the total number of columns. */

double
SparseTensor::getFoldIndexFillFactor() const
{
  int cnt = 0;
  auto start_col = m.begin();
  while (start_col != m.end())
    {
      cnt++;
      const IntSequence &key = start_col->first;
      start_col = m.upper_bound(key);
    }

  return static_cast<double>(cnt)/ncols();
}

/* This returns a ratio of a number of non-zero columns in unfolded
   tensor to the total number of columns. */

double
SparseTensor::getUnfoldIndexFillFactor() const
{
  int cnt = 0;
  auto start_col = m.begin();
  while (start_col != m.end())
    {
      const IntSequence &key = start_col->first;
      cnt += key.getSymmetry().noverseq();
      start_col = m.upper_bound(key);
    }

  return static_cast<double>(cnt)/ncols();
}

/* This prints the fill factor and all items. */

void
SparseTensor::print() const
{
  std::cout << "Fill: "
            << std::fixed << std::setprecision(2) << 100*getFillFactor()
            << std::setprecision(6) << std::defaultfloat << " %\n";
  auto start_col = m.begin();
  while (start_col != m.end())
    {
      const IntSequence &key = start_col->first;
      std::cout << "Column: ";
      key.print();
      auto end_col = m.upper_bound(key);
      int cnt = 1;
      for (auto run = start_col; run != end_col; ++run, cnt++)
        {
          if (cnt % 7 == 0)
            std::cout << "\n";
          std::cout << run->second.first << '(' << run->second.second << ")  ";
        }
      std::cout << "\n";
      start_col = end_col;
    }
}

FSSparseTensor::FSSparseTensor(int d, int nvar, int r)
  : SparseTensor(d, r, FFSTensor::calcMaxOffset(nvar, d)),
    nv(nvar), sym{d}
{
}

void
FSSparseTensor::insert(IntSequence key, int r, double c)
{
  TL_RAISE_IF(!key.isSorted(),
              "Key is not sorted in FSSparseTensor::insert");
  TL_RAISE_IF(key[key.size()-1] >= nv || key[0] < 0,
              "Wrong value of the key in FSSparseTensor::insert");
  SparseTensor::insert(std::move(key), r, c);
}

/* We go through the tensor ‘t’ which is supposed to have single
   column. If the item of ‘t’ is nonzero, we make a key by sorting the
   index, and then we go through all items having the same key (it is its
   column), obtain the row number and the element, and do the
   multiplication.

   The test for non-zero is ‘a != 0.0’, since there will be items which
   are exact zeros.

   I have also tried to make the loop through the sparse tensor outer, and
   find index of tensor ‘t’ within the loop. Surprisingly, it is little
   slower (for monomial tests with probability of zeros equal 0.3). But
   everything depends how filled is the sparse tensor. */

void
FSSparseTensor::multColumnAndAdd(const Tensor &t, Vector &v) const
{
  // check compatibility of input parameters
  TL_RAISE_IF(v.length() != nrows(),
              "Wrong size of output vector in FSSparseTensor::multColumnAndAdd");
  TL_RAISE_IF(t.dimen() != dimen(),
              "Wrong dimension of tensor in FSSparseTensor::multColumnAndAdd");
  TL_RAISE_IF(t.ncols() != 1,
              "The input tensor is not single-column in FSSparseTensor::multColumnAndAdd");

  for (Tensor::index it = t.begin(); it != t.end(); ++it)
    {
      int ind = *it;
      double a = t.get(ind, 0);
      if (a != 0.0)
        {
          IntSequence key(it.getCoor());
          key.sort();

          // check that ‘key’ is within the range
          TL_RAISE_IF(key[0] < 0 || key[key.size()-1] >= nv,
                      "Wrong coordinates of index in FSSparseTensor::multColumnAndAdd");

          auto first_pos = m.lower_bound(key);
          auto last_pos = m.upper_bound(key);
          for (auto cit = first_pos; cit != last_pos; ++cit)
            {
              int r = cit->second.first;
              double c = cit->second.second;
              v[r] += c * a;
            }
        }
    }
}

void
FSSparseTensor::print() const
{
  std::cout << "FS Sparse tensor: dim=" << dim << ", nv=" << nv << ", (" << nr << 'x' << nc << ")\n";
  SparseTensor::print();
}

// GSSparseTensor slicing constructor
/* This is the same as FGSTensor slicing constructor from FSSparseTensor. */
GSSparseTensor::GSSparseTensor(const FSSparseTensor &t, const IntSequence &ss,
                               const IntSequence &coor, TensorDimens td)
  : SparseTensor(td.dimen(), t.nrows(), td.calcFoldMaxOffset()),
    tdims(std::move(td))
{
  // set ‘lb’ and ‘ub’ to lower and upper bounds of slice indices
  /* The same code is present in FGSTensor slicing constructor, see it for
     details. */
  IntSequence s_offsets(ss.size(), 0);
  for (int i = 1; i < ss.size(); i++)
    s_offsets[i] = s_offsets[i-1] + ss[i-1];

  IntSequence lb(coor.size());
  IntSequence ub(coor.size());
  for (int i = 0; i < coor.size(); i++)
    {
      lb[i] = s_offsets[coor[i]];
      ub[i] = s_offsets[coor[i]] + ss[coor[i]] - 1;
    }

  auto lbi = t.getMap().lower_bound(lb);
  auto ubi = t.getMap().upper_bound(ub);
  for (auto run = lbi; run != ubi; ++run)
    {
      if (lb.lessEq(run->first) && run->first.lessEq(ub))
        {
          IntSequence c(run->first);
          c.add(-1, lb);
          insert(c, run->second.first, run->second.second);
        }
    }

}

void
GSSparseTensor::insert(IntSequence s, int r, double c)
{
  TL_RAISE_IF(!s.less(tdims.getNVX()),
              "Wrong coordinates of index in GSSparseTensor::insert");
  SparseTensor::insert(std::move(s), r, c);
}

void
GSSparseTensor::print() const
{
  std::cout << "GS Sparse tensor: (" << nr << 'x' << nc << ")\nSymmetry: ";
  tdims.getSym().print();
  std::cout << "NVS: ";
  tdims.getNVS().print();
  SparseTensor::print();
}
